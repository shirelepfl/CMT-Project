---
title: Lung Cancer Modelisation
author: Shirel Gelfged, Rebecca Levi
geometry: margin=2.5cm
---

The final report should be approximately 5 pages without references and authorship statement, but this is not a hard limit. However, if many figures are generated, they can be included in the appendix to keep the main report centered around the main conclusions. 

With exception to Point 1, the final report should be a standalone document that can be read independently of the project proposal (i.e., your client should be able to understand the project without having to refer back to the original proposal). Therefore, you may reuse (i.e., copy) parts of your project proposal here. The report should be written in a professional manner suitable for submission to a client or publication.

# Deviations from project proposal

In our final project, we made a slight departure from the original proposal while keeping the same overall scope. We modeled two sets of young lungs; one healthy and one representing a smoker; both being invaded by cancer cells. The core of our approach relied on the Fisher-KPP equation. To make the simulations more manageable, we used Neumann boundary conditions and IMEX-ADI schemes, which allowed us to extend our model from 1D to 3D. We did not use the GitHub code provided by Mr. Hernandez. The project followed the timeline outlined in the proposal, and the overall structure and objectives remained consistent with our initial plan.

# Introduction to the problem

Have you ever wondered how quickly cancer cells can spread through the lungs? This question inspired us to build a computational model; written in MATLAB and C; that simulates the 3D invasion of lung tissue by tumor cells. By visualizing this progression, we can illustrate what “early growth” looks like and how rapidly it can advance when nothing slows it down.

We chose to focus on lung cancer because it remains one of the leading causes of cancer-related deaths worldwide. Smoking is the primary risk factor: cigarette smoke contains thousands of chemicals, including many well-known carcinogens that damage the DNA of airway cells. Over time, these mutations accumulate and can trigger cancer. The risk isn’t limited to active smokers, secondhand smoke also affects people nearby, and vaping, although often perceived as safer, still exposes the lungs to toxic compounds.

Raising awareness about lung cancer, especially among younger people, is crucial. Since smoking is a conscious behavior despite widely known consequences, a visual demonstration can have a stronger impact than statistics or warnings alone. Showing how cancer actually spreads inside the lungs may help make the danger more concrete.

To explore this, we modeled two hypothetical 20-year-old men. One leads a healthy lifestyle: regular exercise, good diet, and no smoking. The other is sedentary and smokes heavily, roughly one pack per day (about 20 cigarettes). Comparing these two scenarios allows us to highlight how lifestyle choices can dramatically influence cancer progression in lung tissue.




# Approach used

Describe the approach taken to solve the problem. Include relevant mathematical relationships, models, algorithms, data, etc.  Is the model mechanistic or empirical (e.g., conservation equation, or a parametrized relationship between input and output)? Do you use the program for forecasting/prediction, or inference (e.g., understand model parameters)?

We are using the Fisher-Kolmogorov (Fisher-KPP) equation to model the spatial spreading of cancer cells. This partial differential equation is used to describe how cells evolve over time and space. The equation is written as

$$
\frac{\partial n}{\partial t}
= D\\nabla^2 n + r\\\n\left(1-\frac{n}{K}\right)
$$

where n(x,t) represents the density of cancer cells, D is the diffusion coefficient, r the growth rate, and K the carrying capacity.
The model used is hybrid: the diffusion term is mechanical, representing the physical spread of cells, while the logistic growth term is parametric, capturing empirical characteristics of tumor growth.

To calculate each 
In our model each case starded with the same amount of cancer cells. 
So we began with:
- $r_0 = 0.03 \ \text{day}^{-1}$
- $D_0 = 0.5 \ \text{mm}^2\ \text{day}^{-1}$
- $K_0 = 10^6 \ \text{cells}\ \text{mm}^{-3}$
- $\lambda_0 = 0.01 \ \text{day}^{-1}$

as the initial parameters.
For each case we multiplied these initial parameters by a coefficient corresponding to the change in the evolution of these parameters influenced by their lifestyles.
For the smoker person we used:
- $r_0$ * 1.5 :  smoking and chronic inflammation increase reactive-oxygen species and mutation rate which will increase uncontrolled cell division to a factor ~ 50%; 
- $D_0$ * 1.4 : smoking increases: epithelial-mesenchymal transition (EMT, a process where normal epithelial cells lose their tight connections and become more mobile and invasive), the activity of matrix metalloproteinases (MMP-2 and MMP-9) and as a result, there is more collagen breakdown, which weakens tissue structure making tumor cells more mobile and able to invade lung tissue. A 40% increase matches typical modeling choices in metastasis & EMT studies;
- $K_0$ * 1.3 : Nicotine stimulates angiogenesis, the formation of new blood vessels. As a result, the development of new capillaries is enhanced, which can support tissue repair but also contribute to disease progression. More blood supply induces more oxygen/nutrients which means higher cell density supported. A 30% increase is a literature-aligned estimate; 
- $\lambda_0$ * 0.7 : smoking reduces cytotoxic T cell activity (T cells kill infected or malignant cells presenting specific peptide (MHC I)), NK cell response (Natural Killer (NK) cells are the immune system’s rapid-response bouncers), antigen presentation efficiency (how well a cell chops up proteins into peptides, loads them onto MHC ( I or II), and displays them so T cells can “see” them). Thus, the immune system kills fewer cancer cells. A 30% decrease is biologically plausible.

For the healthy person we used: 
- $r_0$ * 0.8 :  Regular exercise reduces chronic inflammation, improves insulin sensitivity, and enhances apoptosis regulation (the body’s process for eliminating damaged cells), all of which consistently reduce tumor growth signals according to multiple meta-analyses. We use a 20% reduction because of this evidence.
- $D_0$ * 0.9 : In this case, tissues show better collagen organization and less inflammatory remodeling. This means the body produces fewer MMPs (matrix metalloproteinases), enzymes that break down tissue structure. As a result, cell motility (the ability of cells to move) is slightly reduced, helping maintain stable and well-structured tissues. A 10% reduction is appropriate, motility changes are usually subtle.
- $K_0$ * 0.85 : Healthy vasculature and lower angiogenic signaling lead to reduced abnormal blood vessel growth inside tumors, meaning the tumor receives less support. A decrease of about 15% is reasonable and consistent with exercise-related anti-angiogenic effects.
- $\lambda_0$ * 1.4 : Exercise strengthens immune surveillance by increasing NK cell activity, boosting T-cell responses, and improving tumor recognition, which enhances immune-driven tumor killing. A 40% increase matches observed exercise-induced NK cell improvements.
  


For each value, we used sources that are presented in the README file. 

Provide citation to literature where appropriate, particularly to compare your approach to existing work (whether it is similar or different).

# Results

Describe the results. Give your assessment of whether they are reasonable - and how do you determine this?

# Conclusion and outlook

Summarize the approach taken and the answer to the question set out in the problem statement. Describe limitations of the work (outlook) and how it could be improved.

# Authorship statement

Describe contributions of each student on the team to the project in terms of writing of code and report, or project management.

# References

List all references cited in the report in a consistent format. Every reference should have a citation in the text, and every citation in the text should have a corresponding reference.
