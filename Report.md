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

Provide citation to literature where appropriate, particularly to compare your approach to existing work (whether it is similar or different).

# Results

Describe the results. Give your assessment of whether they are reasonable - and how do you determine this?

# Conclusion and outlook

Summarize the approach taken and the answer to the question set out in the problem statement. Describe limitations of the work (outlook) and how it could be improved.

# Authorship statement

Describe contributions of each student on the team to the project in terms of writing of code and report, or project management.

# References

List all references cited in the report in a consistent format. Every reference should have a citation in the text, and every citation in the text should have a corresponding reference.
