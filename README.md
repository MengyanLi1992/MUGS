# MUGS

⚠️ This repository has moved!

The code and R package are now maintained at: https://github.com/celehs/MUGS

We develop MUlti-source Graph Synthesis (MUGS), an algorithm designed to create embeddings for pediatric EHR codes by leveraging graphical information from three distinct sources: (1) pediatric EHR data, (2) EHR data from the general patient population, and (3) existing hierarchical medical ontology knowledge shared across different patient populations. 

Utilizing existing hierarchical medical ontology as prior general knowledge, MUGS facilitates efficient transfer learning by grouping similar codes, thereby enhancing the transferability of knowledge from general to pediatric systems. To address the heterogeneity within code groups and between sites, we propose to decompose a code embedding into three components: the group effect, defined based on the hierarchical medical ontology; the site-nonspecific code effect, capturing characteristics of a code that differ from its group effect and are shared between health systems; and the code-site effect, identifying site-specific traits of a code. Importantly, this decomposition, coupled with penalty functions applied to the code and code-site effects, provides adaptability to varying degrees of heterogeneity within code groups and between sites and protects against the adverse effects of negative knowledge transfer via hyperparameter tuning.

![Flowchart](images/MUGSFlowchart.png)



