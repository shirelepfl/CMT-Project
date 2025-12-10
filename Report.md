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

Why would anyone be interested in this topic: curiosity, societal implications, development of the technology? Provide citation to literature where appropriate.

 Have you ever wondered how quickly cancer cells can spread through the lungs? That question led us to build a computational model, written in MATLAB and C—that simulates how tumor cells invade lung tissue in 3D. By visualizing the spread, we can show what “early growth” looks like and how fast it can advance if nothing slows it down.
 
Why focus on the lungs? Lung cancer is one of the leading causes of cancer death worldwide. Most cases are linked to smoking, because cigarette smoke carries thousands of chemicals—dozens of them are known carcinogens. These substances damage the DNA of cells lining the airways, and over time that damage can build up into cancer. It’s not just cigarettes: secondhand smoke increases risk for people nearby, and vaping is not risk-free—heated aerosols still deliver toxic compounds that can harm the lungs

Definition of project scope:

- what processes will you include / exclude
- what scenarios will you consider

*Note that novelty of the question, approach, or finding is not a requirement for this project.*

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
