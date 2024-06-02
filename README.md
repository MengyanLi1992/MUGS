# MUGS
We develop MUlti-source Graph Synthesis (MUGS), an algorithm designed to create embeddings for pediatric EHR codes by leveraging graphical information from three distinct sources: (1) pediatric EHR data, (2) EHR data from the general patient population, and (3) existing hierarchical medical ontology knowledge shared across different patient populations. 

Utilizing existing hierarchical medical ontology as prior general knowledge, MUGS facilitates efficient transfer learning by grouping similar codes, thereby enhancing the transferability of knowledge from general to pediatric systems. To address the heterogeneity within code groups and between sites, we propose to decompose a code embedding into three components: the group effect, defined based on the hierarchical medical ontology; the site-nonspecific code effect, capturing characteristics of a code that differ from its group effect and are shared between health systems; and the code-site effect, identifying site-specific traits of a code. Importantly, this decomposition, coupled with penalty functions applied to the code and code-site effects, provides adaptability to varying degrees of heterogeneity within code groups and between sites and protects against the adverse effects of negative knowledge transfer via hyperparameter tuning.

![Flowchart](images/MUGSFlowchart.png)

'Main.R' is the R script we used to learn embeddings for MGB and BCH, requiring input of the SPPMI matrices and initial single-site embeddings from these institutions. Since these were constructed based on patient-level data, they could not be made publicly accessible due to privacy constraints. We utilized the existing hierarchical medical ontology of PheCodes, LOINC codes, and RxNorms to group similar codes. For PheCodes, we used the integer level to form groups. For example, PheCode:250.1, PheCode:250.2, and PheCode:250.11 are grouped under PheCode:250. Hierarchical ontologies for LOINC and RxNorm can be found at https://shiny.parse-health.org/hierarchies/. 

After aligning the two sets of embeddings via solving the orthogonal procrustes problem, we trained the initial embeddings for group effects, code effects, and code-site effects by pooling the two sets of embeddings. The process is a standard regression/ANOVE problem, treating code-site effects as residuals. 

We then commenced our core algorithm: updating group effects, code effects, and code-site effects in an alternating and iterative fashion. The three key functions are given in 'MUGSFun.R'. 'GroupEff_par' and 'CodeSiteEff_l2_par'are used to update group effects, and code-site effects, respectively, utilizing parallel computations across multiple cores or machines to enhance speed. 'CodeEff_Matrix' is used to update code effects via matrix computations. 

For hyperparameter tuning, we leveraged silver-standard pediatric PheCode-RxNorm and PheCode-PheCode pairs curated from pediatric articles. This helped select the optimal tuning parameters associated with the penalties on code effects and code-site effects without the need for data splitting. The performance of different sets of embeddings with different tuning parameters was evaluated using 'Embed_Eval_Pediatric.R'. It is designed to assess the accuracy of the embeddings in identifying established related pairs versus random pairs across a wide range of settings.

Although we cannot share the data used to generate the embeddings, we have developed a Shiny App (https://shiny.parse-health.org/multi-view-net/) to support downstream tasks such as feature engineering, construction of pediatric knowledge graphs, and more.



