####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

# Converts Gaussian '09 log file into mol2 file
# - first version June 2012, Andreas Tillack
# - parameters:
#	- Hsummed: if 1, use ESP charges with H summed into connected heavy atom
#	- average_bo:
#		- if 1:
#			- calculate average bond order from same bonds on involved atoms
#			- bond order is rounded to nearest half integer and aromatic type is substituted for b.o.=1.5
#		- if 0 (default):
#			- bond order is calculated from nr of 1/2*(bonding - antibonding electrons of that bond)
#			- bond order is rounded to nearest integer (good for visualization, because it chooses *one* Kekule structure)
#
# CHANGELOG
# - fixed wrong type names for H and F
# - fixed reading of bond order for atoms with more than one letter (i.e. Si)
# - fixed type name for Si (should not have a type, i.e. "Si.2" is wrong)


BEGIN{
# fancy way of initializing a variable with 0 if it does not exist ... (allows to user-define them at execution)
	Hsummed==Hsummed;
	average_bo=average_bo;
# initialize other variables
	ar_spread=0.3;
	nr=0;
	nr_bonds=0;
	charges=0;
	gotESP=0;
	Mulliken=0;
	nbo_count=0;
	found_bonds=0;
}
# Bonds
(($1=="!") && ($2 ~ /^R/)){
	sub(/R/,"",$2);
	bn=$2*1
	sub(/R/,"",$3);
	sub(/\(/,"",$3);
	sub(/\)/,"",$3);
	split($3,b,",");
	bonds[bn][1]=b[1]
	bonds[bn][2]=b[2]
	nr_bonds=bn;
	found_bonds=1;
}
# Bond order from bonding and antibonding orbital occupancies; and lone pairs
((($1=="Natural") && ($2=="Bond")) && ($3=="Orbitals")){
	nbo_count++;
	delete bond_order; # reset bond order array (only use last occurence in log file)
}
# Use bonds from NBO analysis only if there was no bond output from geometry optimization (first in log file)
((nbo_count>=found_bonds) && (($2=="BD") || (($2=="LP") || ($2=="BD*(")))){
	if($2=="BD"){ # bonding orbitals
		A=$6;
		B=$9;
		BO=$10;
		split($4,BC,")");
		if($6=="-"){
			A=$5;
			B=$8;
			BO=$9;
			split($4,BC,")");
		}
		bond_order[A][B]+=BO;
		bond_order[B][A]+=BO;
		if(BC[1]=="1"){
			if(found_bonds==0){
				nr_bonds++;
				bonds[nr_bonds][1]=A
				bonds[nr_bonds][2]=B
			}
			X[A]++;
			X[B]++;
		}
	} else{
		A=$5;
		B=$8;
		C=$6;
		BO=$9;
		split($4,BC,")");
		if($5=="-"){
			A=$4;
			B=$7;
			C=$5;
			BO=$8;
			split($4,BC,")");
		}
		if($2=="BD*("){ # antibonding orbitals
			bond_order[A][B]-=BO;
			bond_order[B][A]-=BO;
		} else{ # lone pair
			if(BC[1]=="1") X[C]++;
		}
	}
}
# Center positions
(($1=="Atomic") && ($2=="Center")){
	result[$3][1]=$6;
	result[$3][2]=$7;
	result[$3][3]=$8;
	if($8==""){ # we have a case like "10.123456-12.123456" or "-2.232323-5.423423" which needs to be fixed
		ncol=split($6,col,"-",seps);
		diff=0;
		if(col[1]==""){ # case like "-2.232323-5.423423"
			diff=1;
		}
		for(i=1+diff; i<=ncol; i++){
			result[$3][i-diff]=seps[i-1]col[i];
#			print $3"."i-diff" = "result[$3][i-diff]
		}
		if(ncol-diff<3){
			ncol2=split($7,col2,"-",seps2);
			for(i=1; i<=ncol2; i++){
				if(seps2[i-1]col2[i]==""){
					ncol--;
				} else{
					result[$3][i+ncol-diff]=seps2[i-1]col2[i];
#					print $3"."i+ncol-diff" = "result[$3][i+ncol-diff]
				}
			}
		}
	}
	nr=$3;
	charges=1;
	if(Hsummed==1) sw=1;
}
# Atom names and charges
((($1=="Mulliken") && ($2=="atomic")) && ($3=="charges:")){
	Mulliken=1;
}
((($1=="Mulliken") && ($2=="charges")) && ($3=="with")){
	if(Hsummed==1) Mulliken=1; else Mulliken=0;
}
((Mulliken==1) && (((NF==3) && (($2+0)!=$2)) && (($2!="Atom") && (gotESP==0)))){ # prefer ESP charges when available, otherwise take Mulliken charges
	atom_charge[$1][1]=$2 # atom name
	atom_charge[$1][2]=$3 # charge
}
(((NF==3) && (($2+0)!=$2)) && ((($1>=1) && ($1<=nr)) && ((charges==1) && ($2!="Atom")))){
	if($1==nr){
		if(sw==0){
			charges=0;
			gotESP=1
		} else sw-=1;
	}
	atom_charge[$1][1]=$2 # atom name
	atom_charge[$1][2]=$3 # charge
}
END{
	if(nr==0){
		print "Gaussian output does not contain necessary information (please use option \"pop=chelp(g)\")."
		exit 1
	}
# Output header
	print "@<TRIPOS>MOLECULE"
	print "*****"
	print " "nr" "nr_bonds" 0 0 0"
	print "SMALL"
	print "USER_CHARGES"
	print ""
	print "@<TRIPOS>ATOM"
# Calculate bond order (if in Gaussian file)
	for(i=1; i<=nr; i++){
		type[i]="";
		if(nbo_count<found_bonds) X[i]=0;
	}
	for(i=1; i<=nr_bonds; i++){
		if(nbo_count>=found_bonds){
			if(bond_order[bonds[i][1]][bonds[i][2]]>0) new_bo[bonds[i][1]][bonds[i][2]]=bond_order[bonds[i][1]][bonds[i][2]]; else new_bo[bonds[i][1]][bonds[i][2]]=2;
		} else{
			new_bo[bonds[i][1]][bonds[i][2]]=2;
			bo[i]=1;
		}
		same_nr[i]=1;
	}
	if(nbo_count>=found_bonds){
		for(i=1; i<=nr_bonds; i++){
			if(average_bo==1){
				for(j=1; j<=nr_bonds; j++){ # Is there a bond to the same atom type?
					same_type=0;
					if(i!=j){
						if(bonds[i][1]==bonds[j][1]){
							if(atom_charge[bonds[i][2]][1]==atom_charge[bonds[j][2]][1]){ # This one is
								same_type=1;
							}
						}
						if(bonds[i][1]==bonds[j][2]){
							if(atom_charge[bonds[i][2]][1]==atom_charge[bonds[j][1]][1]){ # This one is
								same_type=1;
							}
						}
						if(bonds[i][2]==bonds[j][1]){
							if(atom_charge[bonds[i][1]][1]==atom_charge[bonds[j][2]][1]){ # This one is
								same_type=1;
							}
						}
						if(bonds[i][2]==bonds[j][2]){
							if(atom_charge[bonds[i][1]][1]==atom_charge[bonds[j][1]][1]){ # This one is
								same_type=1;
							}
						}
						if(same_type==1){
							new_bo[bonds[i][1]][bonds[i][2]]+=bond_order[bonds[j][1]][bonds[j][2]];
							same_nr[i]++;
						}
					}
				}
			}
			new_bo[bonds[i][1]][bonds[i][2]]/=same_nr[i];
			if(average_bo==1) bo[i]=int(new_bo[bonds[i][1]][bonds[i][2]]+0.5)/2.0; else bo[i]=int(new_bo[bonds[i][1]][bonds[i][2]]/2.0+0.5);
			if(((bo[i]==1.5) && (average_bo==0)) || ((average_bo==1) && (new_bo[bonds[i][1]][bonds[i][2]]/2.0>=(1.5-ar_spread)) && (new_bo[bonds[i][1]][bonds[i][2]]/2.0<=(1.5+ar_spread)))){
				bo[i]="ar";
				type[bonds[i][1]]=".ar";
				type[bonds[i][2]]=".ar";
			}
		}
	}
# Output atom types, positions, and charges
	for(i=1; i<=nr; i++){
		if(((atom_charge[i][1]=="C") || (atom_charge[i][1]=="O")) && ((type[i]!=".ar") && (X[i]>1))){
			if(X[i]<=4){
				if(atom_charge[i][1]=="O"){
					type[i]="."(X[i]);
				} else type[i]="."(X[i]-1);
			} else{
				if(X[i]==5) type[i]=".tbp"; # trigonal bipyramidal
				if(X[i]==6) type[i]=".oh"; # octahedral
			}
		}
		print "\t"i"\t"atom_charge[i][1]"\t"result[i][1]"\t"result[i][2]"\t"result[i][3]"\t"atom_charge[i][1]type[i]"\t1\tLIG1\t"atom_charge[i][2];
	}
# Output connectivity and bond order
	print "@<TRIPOS>BOND"
	for(i=1; i<=nr_bonds; i++){
		print i"\t"bonds[i][1]"\t"bonds[i][2]"\t"bo[i]"\t# "new_bo[bonds[i][1]][bonds[i][2]]/2.0;
	}
}

