---
title: Lung Cancer Modelisation
author: Shirel Gelfged, Rebecca Levi
---

# Deviations from project proposal

In our final project, we made a slight departure from the original proposal while keeping the same overall scope. We modeled two sets of young lungs; one healthy and one representing a smoker; both being invaded by cancer cells. The core of our approach relied on the Fisher-KPP equation. To make the simulations more manageable, we used Neumann boundary conditions and IMEX-ADI schemes, which allowed us to extend our model from 1D to 3D. We did not use the GitHub code provided by Mr. Hernandez. The project followed the timeline outlined in the proposal, and the overall structure and objectives remained consistent with our initial plan.

# Introduction to the problem

Have you ever wondered how quickly cancer cells can spread through the lungs? This question inspired us to build a computational model; written in MATLAB and C; that simulates the 3D invasion of lung tissue by tumor cells. By visualizing this progression, we can illustrate what “early growth” looks like and how rapidly it can advance when nothing slows it down.

We chose to focus on lung cancer because it remains the leading cause of cancer-related deaths worldwide[1]. Smoking is the primary risk factor(85% of all cases [1]): cigarette smoke contains thousands of chemicals, including many well-known carcinogens that damage the DNA of airway cells. Over time, these mutations accumulate and can trigger cancer. The risk isn’t limited to active smokers, secondhand smoke also affects people nearby, and vaping, although often perceived as safer, still exposes the lungs to toxic compounds.

Raising awareness about lung cancer, especially among younger people, is crucial. Since smoking is a conscious behavior despite widely known consequences, a visual demonstration can have a stronger impact than statistics or warnings alone. Showing how cancer actually spreads inside the lungs may help make the danger more concrete.

To explore this, we modeled two hypothetical 20-year-old men. One leads a healthy lifestyle: regular exercise, good diet, and no smoking. The other is sedentary and smokes heavily, roughly one pack per day (about 20 cigarettes). Comparing these two scenarios allows us to highlight how lifestyle choices can dramatically influence cancer progression in lung tissue.




# Approach used

We use the Fisher–Kolmogorov (Fisher-KPP) equation to model the spatial spreading of cancer cells. This partial differential equation describes how cell density evolves over time and space. It is written as:

$$
\frac{\partial n}{\partial t}
= D\\nabla^2 n + r\\\n\left(1-\frac{n}{K}\right)
$$

where n(x,t) represents the density of cancer cells, D is the diffusion coefficient, r the growth rate, and K the carrying capacity.
Our model is hybrid: the diffusion term is mechanical, representing the physical spread of cells, while the logistic growth term is parametric, capturing empirical aspects of tumor development.

In both scenarios, we start with approximately the same initial number of cancer cells, the smoker starts with a little more cancer cells. The baseline parameters we use are:
- $r_0 = 0.03 \ \text{day}^{-1}$
- $D_0 = 0.5 \ \text{mm}^2\ \text{day}^{-1}$
- $K_0 = 10^6 \ \text{cells}\ \text{mm}^{-3}$
- $\lambda_0 = 0.01 \ \text{day}^{-1}$

Each scenario then applies a lifestyle-dependent coefficient to these baseline values, modifying the evolution of the parameters according to the individual’s habits.

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
  
For each coefficient value, we relied on sources listed in the README file.

For our simulation time we chose 4 months so about 120 days; this value corresponds to a good "intermediate growth" time to reach the tumor doubling time [2]. 

The 3D computational grid contains 60 points in each spatial dimension, covering a physical domain of 120 mm in size. The simulation uses a time step of 5 days.

To facilitate our calculations we modeled a Gaussian 3D tumor that starts near the center of our domain, centered in x = L/2. We stated A = $10^5 \ \text{cells}\ \text{mm}^{-3}$ which is a viable maximal density for cells (usually it’s about $10^6$) and sigma = 2 mm, a viable tumor focus. These two values represent well a starting tumor.    

$$
u(x,0) = A\ e^{-\frac{(x - L/2)^2}{2\sigma^2}}
$$

We first apply Neumann boundary conditions to obtain a 1D solution of the KPP equation, imposing zero flux at the domain boundaries ($\frac{\delta_u}{\delta_x}=0$) to ensure that the cancer remains confined within the lungs. We then use the IMEX-ADI method, treating the reaction term implicitly–explicitly in each spatial direction, to extend this solution to 3D. Both methods are described in detail in the comments of our C code.

The MatLab code begins by creating a 3D mesh representing the lungs, modeled as two ellipsoids for the right and left lungs. These are constructed using the meshgrid function for generating coordinates and mathematical equations for the ellipsoids, with the lung volume being defined by boolean masks (mask_right and mask_left) that identify the points within the lung geometry.

Next, the code loads tumor progression data from CSV files containing the number of tumor cells over time for both the healthy individual and the smoker. The readmatrix function is used to extract time steps (t_h, t_s) and tumor cell counts (c_h, c_s) for both cases.

To visualize the tumor growth, the code creates frames at each time step. For each frame, it calculates the number of tumor voxels, in order to make a progression that would be visible yet stay in our domain we had to apply a scale between the size of a 'cell' and the size of our voxels; we used a scale of 14 so 1 voxel is equal to 100 cancer cells. . The save_frame function generates these visualizations by creating a 3D isosurface mesh for the lungs and tumor volume at each time step. The tumor is represented in red, while the lung tissue is drawn with a translucent light pink color. The function also generates and saves each frame as a PNG image and appends these frames to a GIF file that shows the progression of the tumor growth over time.

The script then saves these images as a series of frames in a specified folder and compiles them into a GIF that shows the evolution of cancer in the lungs over time for both the healthy and smoker cases. The imwrite function is used to compile the images into an animated GIF, where each frame represents a snapshot of the tumor's growth.

# Results

<p float="left">
  <img src="https://raw.githubusercontent.com/shirelepfl/CMT-Project/main/images/csv-healthy.png" alt="Healthy lung values" width="20%" />
  <img src="https://raw.githubusercontent.com/shirelepfl/CMT-Project/main/images/csv-smoker.png" alt="Smoker lung values" width="20%" />
</p>

These values, generated by our C code, show how the tumor progression differs between the two subjects. However, this numerical data doesn’t have the same visual impact as the images.

<p float="left">
  <img src="https://raw.githubusercontent.com/shirelepfl/CMT-Project/main/images/Healthy-t0.png"
       alt="Healthy lung at t=0"
       style="width:45%; height:300px; object-fit:contain;" />
  <img src="https://raw.githubusercontent.com/shirelepfl/CMT-Project/main/images/Smoker-t0.png"
       alt="Smoker lung at t=0"
       style="width:45%; height:300px; object-fit:contain;" />
</p>

At the initial time (t=0), both the healthy and smoker lungs have a similar number of tumor voxels, ranging from 155 to 239. The smoker has about 1.5 times more tumor voxels.

<p float="left">
  <img src="https://raw.githubusercontent.com/shirelepfl/CMT-Project/main/images/Healthy-t60.png"
       alt="Healthy lung at t=60"
       style="width:45%; height:300px; object-fit:contain;" />
  <img src="https://raw.githubusercontent.com/shirelepfl/CMT-Project/main/images/Smoker-t60.png"
       alt="Smoker lung at t=60"
       style="width:45%; height:300px; object-fit:contain;" />
</p>

By day 60 (t=60), the difference in tumor progression becomes noticeable. The healthy lung has 1,958 tumor voxels, while the smoker’s lung has 3,935. This shows that the smoker has twice as many tumor cells at this stage.

<p float="left">
  <img src="https://raw.githubusercontent.com/shirelepfl/CMT-Project/main/images/Healthy-t120.png"
       alt="Healthy lung at t=120"
       style="width:45%; height:300px; object-fit:contain;" />
  <img src="https://raw.githubusercontent.com/shirelepfl/CMT-Project/main/images/Smoker-t120.png"
       alt="Smoker lung at t=120"
       style="width:45%; height:300px; object-fit:contain;" />
</p>

At day 120 (t=120), the difference is even more pronounced. The smoker’s lung now has 10,401 tumor voxels, compared to 4,616 in the healthy lung. The smoker’s lung contains 2.25 times as many tumor cells.


From an initial ratio of 1.5x to 2.25x more tumor cells in the smoker, we see a clear, escalating difference. This highlights the far greater progression of tumor growth in the smoker, emphasizing how visualization helps us understand the impact of lifestyle choices on health.

# Conclusion and outlook

In our code, we used the Fisher–Kolmogorov (Fisher-KPP) equation to model the spatial spreading of cancer cells. To solve this equation, we applied Neumann boundary conditions to obtain a 1D solution, which was then extended to 3D using the IMEX-ADI method. The evolution of cancer cells was visualized at each time step in MATLAB as a series of images compiled into GIFs.

These visualizations allow us to clearly demonstrate the impact of lifestyle choices on lung health. From our simulations, it is evident that smoking significantly accelerates cancer cell proliferation: by the end of the simulation, a person who smokes heavily and is in poor physical condition has 2.25 times more cancer cells in the lungs than a healthy, non-smoking individual. This provides a compelling and understandable illustration of the health risks associated with smoking.

Our study does have limitations, however. The mathematical model is complex and relies on several simplifying assumptions, which may affect the precision of the results. It could be improved by using smaller time steps for increased accuracy, incorporating additional biological factors, and validating the simulations against real data. Furthermore, advanced 3D visualizations could make the results more precise and easier to interpret.

# Authorship statement

Our team collaborated closely, meeting for several hours each week. As a result, the project was a collective effort, and the code, text, and other contributions cannot be attributed to any one person.

# References

[1] https://www.who.int/news-room/fact-sheets/detail/lung-cancer
[2] https://pubmed.ncbi.nlm.nih.gov/8072198/ 
The rest of the citations that we mentioned for the approximation made can be found in the README file. 
