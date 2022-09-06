#for delta in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1;do
#for i in `seq 0 999`; do
#echo $delta, $i
#cp /home/mrivas/jupyter/projectPheFlux/results/validation/simulatedResults/predFluxes/Flux_Ecoli_noATPM_fixedBiomass_noCap_iJO1366_predictedPheflux_shuffled_False_delta_${delta}_coverage_1.0_condition_${i}.txt Ecoli_iJO1366_predictedPheflux_Lambda_${delta}_replicate_${i}.txt
#done; done

# Copie FBAl2 predictions

for delta in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1;do
#cp /home/mrivas/jupyter/projectPheFlux/results/validation/simulatedResults/predFluxes/Flux_Ecoli_noATPM_fixedBiomass_noCap_iJO1366_predictedFBAl2_shuffled_False_delta_${delta}_coverage_1.0_condition_0.txt Ecoli_iJO1366_predictedFBAl2_Lambda_${delta}.txt
#mv coli_iJO1366_predictedFBAminl2_Lambda_${delta}.txt Ecoli_iJO1366_predictedFBAminl2_Lambda_${delta}.txt
done
