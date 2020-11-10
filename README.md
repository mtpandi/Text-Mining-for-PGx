# Text Mining the PGx literature

### Project description
Text-mining in biomedical literature is an emerging field which has already shown to have a variety of implementations in many research areas, including genetics, personalized medicine, and pharmacogenomics. In this study, we describe a novel text-mining approach for the extraction of pharmacogenomics associations. The code that was used towards this end was implemented using R programming language, either through custom scripts, where needed, or through utilizing functions from existing libraries. Articles (abstracts or full texts) that correspond to a specified query were extracted from PubMed, while concept annotations were derived by PubTator Central. Terms that denote a Mutation or a Gene as well as Chemical compound terms corresponding to drug compounds were normalized and the sentences containing the aforementioned terms were filtered and preprocessed to create appropriate training sets. Finally, after training and adequate hyperparameter tuning, four text classifiers were created and evaluated (FastText, Linear kernel SVMs, XGBoost, Lasso and Elastic-Net Regularized Generalized Linear Models) with regards to their performance in identifying pharmacogenomics associations. Although further improvements are essential towards a proper implementation of this text-mining approach in the clinical practice, our study stands as a comprehensive, simplified and up-to-date approach for the identification and assessment of research articles enriched in clinically relevant pharmacogenomics relationships. Furthermore, this work highlights a series of challenges concerning the effective application of text mining in biomedical literature, whose resolution could substantially contribute to the further development of this field.  


![Alt text](https://github.com/mtpandi/Text-Mining-for-PGx/blob/master/PROJECTlayout.png)


### References
Pandi M-T, van der Spek PJ, Koromina M and Patrinos GP (2020) A Novel Text-Mining Approach for Retrieving Pharmacogenomics Associations From the Literature. Front. Pharmacol. 11:602030. doi: 10.3389/fphar.2020.602030
