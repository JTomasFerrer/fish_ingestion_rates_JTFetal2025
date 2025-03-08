This reprository contains all the code related to our work "Ultrasonography as a non-invasive 
technique for assessing diet effects on fish gonad dynamics", by Joaquim Tomàs-Ferrer, Irene 
Moro-Martínez, Enrique Massutí-Pascual, Amàlia Grau and Miquel Palmer.

The R script jtf_weight_dynamics (and the input data file input.RData) are designed for analysing the effect of diet on fish gonad 
dynamics using ultrasonography data. This script contains several steps:

1. Loading data and libraries

2. Stan model specification
The script writes a Stan model into a file named model.stan, this model is structured to estimate
the relationship between food intake and 3 key parameters related to gonad dynamics:
-mu: the time at which maximum gonad weight is attained.
-sigma: the spread of the reproductive period.
-h: the maximum gonad weight.
The model includes both fixed effects (overall relationships with food intake) and random effects
(variability between tanks and fish).

3. Compiling the model

4. Initialising the model

5. Running the model

6. Saving the results

7. Results summary
The main purpose of this code is to estimate how diet influences the timing, duration and
magnitude of gonad development in fish, as observed through ultrasonography. The Bayesian model
accounts for both individual fish differences and tank-level differences, allowing a robust
analysis of the data. The results can help understanding the relationship between food intake
and reproductive traits, which is important for both fisheries management and aquaculture.
