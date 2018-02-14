####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

BEGIN{
	items=0
	start_column=3
	end_column=start_column
	zeroed=0
}
(!/^#/){
	end_column=NF
	if(zeroed==0){
		for(i=start_column; i<=NF; i++){
			sum[i] = 0
			sum2[i] = 0
		}
		zeroed=1
	}
	for(i=1; i<start_column; i++){
		description[i]=$i
	} 
	for(i=start_column; i<=NF; i++){
		sum[i] = sum[i]+$i
		sum2[i] = sum2[i]+$i*$i
	}
	items++
}
END{
	for(i=1; i<start_column; i++) printf("%s\t", description[i])
	for(i=start_column; i<=end_column; i++){
		variance=sum2[i]/items-(sum[i]/items)**2;
		if(variance<1E-10) variance=0.0
		printf("%f\t%f", sum[i]/items,sqrt(variance))
		if(i<end_column) printf("\t")
	}
}
