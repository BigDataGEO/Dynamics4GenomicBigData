% This function receives a row vector ([mean_curve]) of T elements, representing a gene's gene
% expression over T time points, and returns a matrix of size [number_of_curves]xT where each row
% is a normal deviate of [mean_curve] with standard deviation [standard_deviation].

function simulated_gene_expression = generate_simulated_curves_from_mean_curve(number_of_curves, mean_curve, standard_deviation)
  
  simulated_gene_expression = [];

  for i=1:number_of_curves
    simulated_gene_expression = [simulated_gene_expression; normrnd(mean_curve, standard_deviation)];
  end
end