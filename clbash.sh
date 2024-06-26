#! usr/bin/bash

#BASH RUNS SCRIPT FOR COSMOLATTICE
#This script allows for repeated runs of CosmoLattice and a simple scanning of the parameter space.
#In this specific case, the parameters to be scanned are lambda and nu, while the coupling g is fixed.
#It's a bit rough but it gets the job done.


# High Intensity Colours
IBlack='\033[0;90m'       # Black
IRed='\033[0;91m'         # Red
IGreen='\033[0;92m'       # Green
IYellow='\033[0;93m'      # Yellow
IBlue='\033[0;94m'        # Blue
IPurple='\033[0;95m'      # Purple
ICyan='\033[0;96m'        # Cyan
IWhite='\033[0;97m'       # White

Colour_Off='\033[0m'      # Text Reset


printf "\n"
printf "\n"
echo -e "${ICyan}"
echo -e "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
echo -e "XXXX                                             XXXX"
echo -e "XXXX        COSMOLATTICE BASH RUNS SCRIPT        XXXX"
echo -e "XXXX                                             XXXX"
echo -e "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX${Colour_Off}"
printf "\n"
printf "\n"

# DEFAULT MODEL PARAMETERS

Hkin_val=10^5
# Self-coupling parameter
log10_lambda_min=4.5
log10_lambda_max=5.2

#xi_min=1
#xi_max=50

lambda_min=$(echo "scale=7; e($log10_lambda_min*l(10))" | bc -l)
lambda_max=$(echo "scale=7; e($log10_lambda_max*l(10))" | bc -l)

# Scanning steps of the self-coupling parameter
lambda_steps=8
#xi_steps=5

log10_m_min=4
log10_m_max=7

m_min=$(echo "scale=7; e($log10_m_min*l(10))" | bc -l)
m_max=$(echo "scale=7; e($log10_m_max*l(10))" | bc -l)


# Scanning steps of the self-coupling parameter
m_steps=5

# Lattice parameters 
n_lattice=128
t_max=2000


#INTERFACE AND PARAMETERS

echo -e "${ICyan}"
echo -e "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
echo -e "XXXX                 RUNS SETTINGS               XXXX"
echo -e "XXXX                                             XXXX"
echo -e "XXXX DEFAULT MODEL PARAMETERS:                   XXXX"
echo -e "XXXX                                             XXXX"
echo -e "XXXX  1) log10(lambda min) = $log10_lambda_min         "
echo -e "XXXX  2) log10(lambda max) = $log10_lambda_max         "
echo -e "XXXX                                                 "
echo -e "XXXX  3) lambda min = $lambda_min                       "
echo -e "XXXX  4) lambda max = $lambda_max                       "
echo -e "XXXX                                                 "
echo -e "XXXX  5) lambda steps = $lambda_steps                   "
echo -e "XXXX                                                 "
echo -e "XXXX  6) m min = $m_min                               "
echo -e "XXXX  7) m max = $m_max                               "
echo -e "XXXX                                                 "
echo -e "XXXX  8) m steps = $m_steps                           "
echo -e "XXXX                                                 "
echo -e "XXXX DEFAULT LATTICE PARAMETERS:                 XXXX"
echo -e "XXXX                                             XXXX"
echo -e "XXXX  9) N LATTICE = $n_lattice                                 "
echo -e "XXXX  10) MAXIMUM TIME = $t_max                          "
echo -e "XXXX                                             XXXX"
echo -e "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
echo -e "${Colour_Off}"


# Interactive change of default parameters

printf "\n"
echo -e "Do you want to proceed with default settings? (y/n)"
read def_set

answ="y"

if [ "$def_set" == "n" ];
then
	while [ "$answ" == "y" ]
	do
	printf "\n"
	echo -e "Which parameter would you like to change? Insert number"
	read param_num
	
	case "$param_num" in 
		1)
		echo -e "1) log10(lambda min) = "
		read new_param
		log10_lambda_min=$new_param
		;;

		2)
		echo -e "2) log10(lambda max) = "
		read new_param
		log10_lambda_max=$new_param
		;;

		3)
		echo -e "3) lambda min = "
		read new_param
		lambda_min=$new_param
		;;

		4)
		echo -e "4) lambda max = "
		read new_param
		lambda_max=$new_param
		;;

		5)
		echo -e "5) lambda steps = "
		read new_param
		lambda_steps=$new_param
		;;

		6)
		echo -e "6) m min = "
		read new_param
		m_min=$new_param
		;;

		7)
		echo -e "7) m max = "
		read new_param
		m_max=$new_param
		;;

		8)
		echo -e "8) m steps = "
		read new_param
		m_steps=$new_param
		;;

		9)
		echo -e "10) N LATTICE = "
		read new_param
		n_lattice=$new_param
		;;

		10)
		echo -e "11) MAXIMUM TIME = "
		read new_param
		t_max=$new_param
		;;


		*)
		echo -e "Unknown"
		;;

	esac
	
	printf "\n"
	echo -e "Want to change something else? (y/n)"
	read more_changes
	answ=$more_changes
	
	done
	
	printf "\n"
	echo -e "${ICyan}"
	echo -e "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
	echo -e "XXXX                 RUNS SETTINGS               XXXX"
	echo -e "XXXX                                             XXXX"
	echo -e "XXXX DEFAULT MODEL PARAMETERS:                   XXXX"
	echo -e "XXXX                                             XXXX"
	echo -e "XXXX  1) log10(lambda min) = $log10_lambda_min         "
	echo -e "XXXX  2) log10(lambda max) = $log10_lambda_max         "
	echo -e "XXXX                                                 "
	echo -e "XXXX  3) lambda min = $lambda_min                       "
	echo -e "XXXX  4) lambda max = $lambda_max                       "
	echo -e "XXXX                                                 "
	echo -e "XXXX  5) lambda steps = $lambda_steps                   "
	echo -e "XXXX                                                 "
	echo -e "XXXX  6) m min = $m_min                               "
	echo -e "XXXX  7) m max = $m_max                               "
	echo -e "XXXX                                                 "
	echo -e "XXXX  8) m steps = $m_steps                           "
	echo -e "XXXX                                                 "
	echo -e "XXXX DEFAULT LATTICE PARAMETERS:                 XXXX"
	echo -e "XXXX                                             XXXX"
	echo -e "XXXX  9) N LATTICE = $n_lattice                                 "
	echo -e "XXXX  10) MAXIMUM TIME = $t_max                          "
	echo -e "XXXX                                             XXXX"
	echo -e "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
	echo -e "${Colour_Off}"
	printf "\n"

fi

# Selection of the output folder
echo -e "Insert name of the folder (write nnumber for nmc_runs_n):"
read foldername

cd /home/tomas_mendes/projects/cosmolattice
mkdir -p nmc_runs_$foldername
cd nmc_runs_$foldername
make clean-cmake
cmake -DMODEL=GB_wMass -DMPI=ON ../
make cosmolattice

# Writing the file containing the lambda and nu values
#echo "lambda       nu" >> parameter_scan.txt

# LOOPS FOR PARAMETER SCANNING
for (( i=0; i<=lambda_steps; i++ )) do #loop over lambda values
	
	#Current value of lambda in the run
	lambda_exp=$(echo "scale=10; $log10_lambda_max + $i*($log10_lambda_min - $log10_lambda_max)/($lambda_steps)" | bc -l)
	lambda_val=$(echo "scale=10; e($lambda_exp*l(10))" | bc -l)

	#xi_val=$(echo "scale=3; $xi_min + $i*($xi_max-$xi_min)/$xi_steps" | bc )

	if [ $lambda_steps == 0 ];
		then
		lambda_val=0
	fi

	for (( j=0; j<=m_steps; j++ )) do #loop over nu values
	
		# Current value of nu in the run
		m_exp=$(echo "scale=10; $log10_m_max + $j*($log10_m_min - $log10_m_max)/($m_steps)" | bc -l)
		m_val=$(echo "scale=10; e($m_exp*l(10))" | bc -l)
		#m_val=$(echo "scale=3; $m_min + $j*($m_max-$m_min)/$m_steps" | bc )

		# Minimum and maximum tachyonic momenta
		k_min=$(echo "scale=7;  $lambda_val / ($Hkin_val*sqrt(48))" | bc )
		#k_min=$(echo "scale=7;  1 / (sqrt(48*$xi_val))" | bc )
		k_max=$(echo "scale=7; 1" | bc )

		# Run counter
		n_run=$(echo "$i * ($m_steps + 1) + 1 + $j" | bc -l)
		tot_run=$(echo "scale=3; ( 1 + $lambda_steps ) * ( 1 + $m_steps )" | bc )

		# Creation of the numbered output folders 
		printf -v dir_index "%05d" $n_run
		mkdir -p run_$dir_index

		# Writing the current parameter combination into a text file
		echo "$lambda_val $m_val $dir_index" >> parameter_scan.txt

		deltat=0.1
		if [ $k_min < 0.1 ];
			then
			deltat=0.01
		fi

		# PARAMETERS IN THE CURRENT RUN

		printf "\n"
		printf "\n"
		echo -e "${IGreen}"
		echo -e "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
		echo -e "XXXX                                                 "
		echo -e "XXXX         COSMOLATTICE RUNNING                    "
 		echo -e "XXXX           RUN $n_run OF $tot_run                "
		echo -e "XXXX                                                 "
		echo -e "XXXX CURRENT PARAMETER VALUES IN THE RUN:            "
		echo -e "XXXX                                                 "
		echo -e "XXXX   log10(lambda) = $lambda_exp                   "
		echo -e "XXXX          lambda = $lambda_val                   "
		echo -e "XXXX   log10(M) = $m_exp                   "
		echo -e "XXXX              M = $m_val                         "
		echo -e "XXXX												  "
		echo -e "XXXX			k_min = $k_min						  "
		echo -e "XXXX			k_max = $k_max						  "
		echo -e "XXXX                                                 "
		echo -e "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
		echo -e "${Colour_Off}"
		printf "\n"
		printf "\n"

		mpirun GB_wMass input=../src/models/parameter-files/GB_wMass.in N=$n_lattice kIR=$k_min dt=$deltat tMax=$t_max tOutputVerb=10 kCutOff=$k_max M=$m_val cutoff=$lambda_val baseSeed=1922 outputfile=./run_$dir_index/

		# Launching CosmoLattice
        #sbatch --dependency=singleton --job-name=GB_DM /home/tomas_mendes/projects/cosmolattice/cl_GB_bash.sh $n_lattice $t_max $Hkin_val $cutoff $m_val  $k_min $k_max $n_run $dir_index


	done
done
