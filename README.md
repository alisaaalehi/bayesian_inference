# bayesian_inference

Bayesian procedure for regression based on latent variable approach.
In this project, I am going to implement Gibbs sampling to estimate parameters of a regression problem with specific priors.
Let's look at the problem itself:



We have Dose-Response Data from an Insecticide Study in which insects were exposed to various dose levels of an insecticide and the following is the result of this:


Dose	 N	 A
1.69	 60	6
1.72	 62	13
1.76	 63	20
1.78	 60	30
1.81	 64	53
1.84	 60	55
1.86	 62	61
1.88	 64	62


where N = number exposed; A = number adversely affected.
Let P(d) = P( Exposed insect to dose d is adversely affected). and let P(d) =  Φ(β0 + β1 x d + β2 x d^2; 1),
the goal is finding the Bayesian estimates of β0, β1, β2 under certain assumptions explained in Project.pdf file.




- Implementation is in R
- For a detailed explanation of the project definition please look at the Project.pdf file
- the file Report.pdf contains results and mathematical derivations of the solutions
- This is part of a project for Bayesian Inference class under Dr. E. Olúṣẹ́gun George.
- The implementations are designed based on tutorials by Andreas C. Kapourani (https://rstudio-pubs-static.s3.amazonaws.com/208180_b659633007eb45aa9c48e4c50b8afc07.html)
