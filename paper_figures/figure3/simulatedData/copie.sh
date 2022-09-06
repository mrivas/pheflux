
origin="/home/mrivas/jupyter/projectPheFlux/data/validation/simulatedData/"
destiny="/home/mrivas/jupyter/projectPheFlux/paper/paper_figures/figure3/simulatedData/"

for delta in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do
cp ${origin}Ecoli_noATPM_fixedBiomass_iJO1366_simulatedFluxes_delta_${delta}.txt ${destiny}Ecoli_iJO1366_simulatedFluxome_Lambda_${delta}.txt
cp ${origin}Ecoli_noATPM_fixedBiomass_iJO1366_simulatedFPKM_delta_${delta}.txt ${destiny}Ecoli_iJO1366_simulatedTranscriptome_Lambda_${delta}.txt
done
