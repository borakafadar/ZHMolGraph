communications biol ogy 

Article 

A Nature Portfolio journal 

https://doi.org/10.1038/s42003-025-07694-9 

## RNA-protein interaction prediction using network-guided deep learning 

Check for updates 

Haoquan Liu[1] , Yiren Jian © 2, Chen Zeng © 3 & Yunjie Zhao © 1 & 

Accurate computational determination of RNA-protein interactions remains challenging, particularly when encountering unknown RNAs and proteins. The limited number of RNAs and their flexibility constrained the effectiveness of the deep-learning models for RNA-protein interaction prediction. Here, we introduce ZHMolGraph, which integrates graph neural network and unsupervised large language models to predict RNA-protein interaction. We validate ZHMolGraph predictions on two benchmark datasets and outperform the current best methods. For the dataset of entirely unknown RNAs and proteins, ZHMolGraph shows an improvement in achieving high AUROC of 79.8% and AUPRC of 82.0%. This represents a substantial improvement of 7.1%–28.7% in AUROC and 4.6%–30.0% in AUPRC over other methods. We utilize ZHMolGraph to enhance the challenging SARS-CoV-2 RPI and unbound RNA-protein complex predictions. Such enhancements make ZHMolGraph a reliable option for genome-wide RNA-protein prediction. ZHMolGraph holds broad potential for modeling and designing RNA-protein complexes. 

RNA–proteincomplexesare essentialfor manycellularprocesses,including gene transcription and post-transcriptional gene regulation[1][–][3] . It is also worth noting that RNA–protein interactions are closely related to diseases[4] . For example, many retroviruses, like HIV, heavily rely on RNA–protein interactions to facilitate their replications within humans[5][,][6] . As RNA – mutates rapidly, determining the structure of RNA protein complexes in time can be challenging. Therefore, an urgent need is to rapidly identify RNA’s protein binding partners leveraging the existing RNA–protein interaction networks and sequencing data. 

Experimental techniques such as X-ray crystallography, nuclear magnetic resonance (NMR), and cryo-electron microscopy are used to determine the structure of RNA–protein complexes and their interactions[7][–][10] . High-throughput techniques, including PAR-CLIP, RNAcompete, RIP-Chip, and HITS-CLIP, are used to measure RNA–protein interactions without requiring structures[11][–][14] . However, they all suffer the drawback of being time-consuming and costly. Some computational methods can predict the RNA–protein complex structures through homologous fragment modeling, docking, or molecular dynamics simulation with the help of the structure scoring functions[15][–][18] . Unfortunately, due to the limited experimental structures available for RNA and protein[19] , the current scoring function needs more improvements to determine the binding conformation of recently discovered RNA–protein complexes accurately. Alphafold3 can be applied to predict the RNA–protein structures[20] . However, it is more effective and accurate to predict the existence of interactions between RNA and protein before the complex 

structuralmodeling.Consequently,identifyingpairsofinteractionsbetween RNA and protein at a genome-wide level is a crucial step in understanding cell regulatory mechanisms. 

State-of-the-art RNA–protein interaction (RPI) prediction methods primarily rely on traditional machine learning and deep learning techniques. One such approach is RPIseq, which encodes 4-mer RNA and 3-mer protein sequences as feature vectors for input into support vector machine (SVM) and random forest (RF) classifiers. These classifiers were employed to predict the presence or absence of RNA–protein interacting pairs[21] . Traditional machine learning methodologies depend on established techniques,suchasRFandSVM,whichexhibiteffectivetrainingcapabilitiesand computational efficiency. Nevertheless, considering the exponential expansion of RNA and protein sequence data, these methodologies have proven insufficient. Furthermore, deep learning has contributed to advan– cing RNA protein interaction (RPI) prediction. 

Deep learning models are used to learn the patterns from the node degrees in the RPI network. For example, IPMiner has started to contribute to RPI prediction by using stacked autoencoders to extract latent features of high-level abstraction from RNA and protein K-mer sequence vectors, thus uncovering more features within the sequences[22] . Graph neural networks have also been employed to integrate information within nodes to predict RPIs based on the available RPI network. The NPI-GNN framework combines Graph neural network and top-k pooling within the SEAL framework to address NPI prediction tasks[23][,][24] . SEAL, a graph neural network-based framework, constructs subgraphs from links to effectively reframe the link 

> 1Institute of Biophysics and Department of Physics, Central China Normal University, Wuhan, 430079, China. 2 College, Hanover, NH, 03755, USA.[3] Department of Physics, The George Washington University, Washington, DC, 20052, USA. 

&) e-mail: yjzhaowh@ccnu.edu.cn 

Communications Biology | (2025) 8:247 

1 

Article 

https://doi.org/10.1038/s42003-025-07694-9 

prediction challenge into a binary classification problem based on subgraph characteristics. However, some imbalances exist in the distribution of connected neighbors associated with different RNAs and proteins. The feature vector of nodes also causes an imbalance in the deep learning models. They fail to learn the binding propensity of “orphan” RNAs and proteins with few or no connections or homologs in the dataset[25][,][26] . Deep learning models capture annotation imbalances, leading to imbalanced predictions of RNA–protein binding. Therefore, the accuracy of deep learning methods in predicting interactions between unknown RNA and proteins is low. 

In our research, we strive to enhance the accuracy of predicting RNA–protein interactions, especially when it comesto unknown RNAs and proteins. Our approach characterizes the RPI networks constructed from various sources at the entire biomolecule (protein or RNA) and finer residue/nucleotide scales. We then introduce ZHMolGraph, which combines the graph neural network with unsupervised large language models (LLMs) to address the annotation imbalance inherent in the existing RPI networks. - Compared to traditional machine learning techniques, ZHMolGraph effi ciently integrates information from known network interactions to enhance accuracy in benchmark testing. ZHMolGraph utilizes unsupervised learning to effectively capture insights from the RNA-FM and ProtTrans LLMs[27][,][28] . With the help of those LLMs on sequencing-scale data, ZHMolGraph overcomes limitations posed by limited binding data, enhancing its generalizability to unknown and protein pairs. The improvement of the ZHMolGraph model exhibits potential in aiding the – prediction of RNA protein structures. 

## Results 

## The overview 

– This work aims to improve RNA Protein RPI prediction. In the first step, we constructed three distinct networks using structural, high-throughput, 

and literature-mining validated data to understand the characteristics of RPIs (see Fig. 1). We break down the interactions in a structural RNA–protein complex to those between their residues and nucleotides. A pair of residue and nucleotide is an interaction if they exhibit non-covalent interactions within a cutoff distance of 8 Å. We also extract RNA–protein interactions from the RNAInter database built from high-throughput techniques and the NPInter5 database built from literature mining. See the “Methods” section for more detailed information on constructing the RPI networks. Our analysis of RPI networks from multiple sources has established common characteristics in biomolecule and finer residue-level RPI networks. 

In the second step, we propose a deep-learning pipeline called ZHMolGraph that utilizes a network-sampling learning strategy with unsupervised LLMsnodefeaturestooptimizetheexplorationofthebinding properties of RNAs and proteins. The workflow of ZHMolGraph involves generating embedding features for RNA and protein sequences using RNAFM and ProtTrans (see Fig. 2). The LLMs embedding features are subsequently fed into the graph neural network model to integrate and aggregate the network information from the RPI network. Lastly, the LLMs embedding and graph neural network sampling features are concatenated and input into the VecNN to predict the binding likelihood. 

## The characteristics of RPI networks 

We have established three RPI networks. The first is the structural network (seeFig. 1aandSupplementaryData1),whichcontains1198RNAand3399 protein nodes connected by 7699 edges representing their interactions. The second is the high-throughput interactions network (see Fig. 1b and Supplementary Data 2), which consists of 23,205 RNA and 4713 protein nodes linked by 82,170 edges. Finally, the literature-mining validated network (see Fig. 1c and Supplementary Data 3) has 1937 RNA and 1716 protein nodes 

Fig. 1 | A workflow to construct RPI networks. a The structural RNA–protein interaction network comprises 1198 RNAs and 3399 proteins from the PDB database. b The high-throughput sequencing RNA–protein interaction network encompasses 23,205 RNAs and 4713 proteins from the RNAInter database. c The literature-mining validated RNA–protein interaction network consists of 1937 RNAs and 1716 proteins retrieved from the NPInter database. In (a–c), RNAs and proteins are denoted by pink and green nodes. The solid links indicate interactions between RNAs and proteins. d A – selected subgraph from the structural RNA protein interaction network. e The RNA–protein residue network is derived from the subgraph (d), and it includes 2125 nucleotides and 3907 amino acid residues. 

Communications Biology | (2025) 8:247 

2 

Article 

https://doi.org/10.1038/s42003-025-07694-9 

Fig. 2 | Overview of ZHMolGraph framework. a ZHMolGraph pipeline leverages unsupervised LLMs methods to generate embeddings for RNAs and proteins. These embeddings are then passed on to a graph neural network, which produces network embeddings. Both the LLMs and network embeddings are used to train deep models. b ZHMolGraph architecture uses RNA-FM and ProtTrans to generate the LLMs embeddings. It also utilizes the graph neural network algorithm to sample the network embedding of nodes. On benchmark datasets, the VecNN is trained in a 5-fold cross-validation setup. The neural network’s final output is averaged prediction over the five folds. 

connected by 4329 edges. These networks’ node and edge counts are crucial for gaining insights into the complexities of RPI networks. 

## Scale-free with high modularity in the molecule RPI networks 

We analyzed the topology of structural networks, validated networks based on literature mining, and high-throughput networks. We found a fat-tailed distribution of interactions associated with all nodes, proteins, and RNAs in the structuralnetwork(see Fig. 3a,c,e). The degreedistributionfor allnodes in the structural network follows a power-law characterized by a degree exponent (γ) of 2.561 (see Fig. 3b, d, and f on the double logarithmic axes). The exponent for node degrees is 2.135 for RNA and 3.203 for protein. This – suggests that the RNA protein structural interaction network is scale-free, withmostproteinsandRNAsshowingonlylimitedinteractions.However,a select few hub nodes possess an exceptionally high number of binding records. We also found a similar fat-tailed distribution in the highthroughput and literature-mining validated networks. The degree exponents for the high-throughput networks are 2.665, 2.361, and 3.139 for all nodes, RNA nodes, and protein nodes, respectively (see Supplementary – Fig.1a f).Thedegreeexponentsintheliterature-miningvalidatednetworks are 2.537, 2.629, and 2.543 for all nodes, RNA nodes, and protein nodes, – respectively (see Supplementary Fig. 2a f). 

The study shows that RPI networks have a scale-free topology. It highlights differences in the connectivity preferences of nodes within the network. Hubs are the nodes that have connections with many other nodes, whilemost nodesconnectonlywith afewothers.The largelanguagemodels interpret an entire sequence as one sentence and its amino acids (nucleotides) as individual words. The connectivity pattern suggests that the formation of RPIs is not random but follows specific rules. Graph neural networks can potentially improve the “sentence-sentence” attention mechanism. 

The topological coefficient measures how extensively a node in a network shares neighbors with other nodes. We calculate this coefficient for nodes in the network that have common neighbors with other nodes. We also investigate the correlation between the topological coefficient and node degree,showinganapparentpower-law decay.Inthestructuralnetwork,we observed an anti-correlation between the degree of all nodes and the topological coefficient. Supplementary Figs. 3a, c, and e depict this 

relationship with a Spearman correlation coefficient of all nodes, RNA nodes, and protein nodes are −0.927, −0.856, and −0.944, respectively. A power-law decay is also identified in Supplementary Figs. 3b, d, f. 

High-throughput and literature-mining validated networks also exhibit an anti-correlated and power-law decay pattern. A significant anti− correlatedSpearmancorrelationof 0.977wasobservedbetweenthedegree and topological coefficient of all nodes in the high-throughput networks. In addition, RNA and protein nodes exhibit anti-correlated Spearman correlation coefficients of −0.998 and −0.841, respectively (see Supplementary – Fig. 4a f). Similar results were observed in the literature-mining validated networks. For all nodes, RNA nodes, and protein nodes, the Spearman correlation coefficients between the degree and topological coefficient are −0.982, −0.985, and −0.975, respectively (see Supplementary Fig. 5a–f). These findings suggest a strong relationship between network nodes’ degree and topological coefficient. The nodes with high degrees generally do not have a significantly greater number of common neighbors than nodes with lower degrees across all three types of networks, suggesting a high level of modularity in these networks. 

After analyzing the network properties of various sources, it has been observed that the RPI network displays scale-free characteristics and high modularity. These findings suggest a sense of order and consistency within the complex RPI network and provide a fresh perspective on its intricate topological structure. RPI complex network systems display well-defined connectivity relationships, allowing large language models to map associations between phrases and infer interactions between complexes. 

## Scale-free with high modularity in the residue RPI network 

– To investigate the residue-level RNA protein network, we focus on a sub– graph of the SARS-CoV-2-Nsp1-40S complex in the RNA protein structural interaction network (PDB code: 6ZOL), including an RNA node with 58 protein nodes connected by 58 edges. Thus, we generate a residue-level network with 2125 RNA nucleotides and 3907 protein residues (see Fig. 1d, e and Supplementary Data 4). 

We analyzed the distribution of degrees in the residue network and found a fat-tail distribution for all nodes. This indicates the presence of a scale-free phenomenon in RPIs at the residue level (see Supplementary Figs. 6a, c, e). The degree exponents for all nodes, RNA nucleotide nodes, 

Communications Biology | (2025) 8:247 

3 

Article 

https://doi.org/10.1038/s42003-025-07694-9 

– Fig. 3 | Degree distribution of RNA protein structural interaction network from PDB database. a, c, e The degree distributions of a all nodes, c RNA nodes, and e protein nodes in RNA–protein structural interaction network. b, d, f The degree distributions of b all nodes, d RNA nodes, and f protein nodes in the RNA–protein structural interaction network are shown in double logarith– mic axes (log log plot), indicating that the degree probability of all nodes, RNA nodes, and protein nodes are well approximated by a power law with degree exponents 2.561, 2.135, and 3.203, respectively. Source data are provided in Supplementary Data 6. 

and protein residue nodes are 2.790, 2.838, and 4.876, respectively (see Supplementary Fig. 6b, d, f). We also analyzed the correlation between the degree and topological coefficient of nodes within the residue network. The results show an anti-correlation between degree and topological coefficient for all nodes. The Spearman coefficients for the degree and topological coefficients of all nodes, RNA nodes, and protein nodes are −0.938, −0.929, and −0.952, respectively (see Supplementary – Fig. 7a f). 

Looking at the RPI network at the residue level offers a more detailed and comprehensive understanding. Surprisingly, our observations at this smaller scale align with those at the biomolecular level, showing that the network has scale-free characteristics and high modularity. This reveals the existence of binding hotspots within the sequence and suggests that the topological features within the RPI network are consistent, even at a finer residue level. This residue-level information holds the potential to enhance the prediction of RPI. 

## Attachment ability of nodes in the structural RPI network 

Our study on the interaction networks between RNA and protein molecules reveals that these networks exhibit scale-free and high modularity features. This means that RPI attaches preferentially to a limited number of nodes with a substantial number of connections. In contrast, the majority of nodes only have a few connections. Therefore, we aim to investigate whether the 

RNA–proteinnetworkexhibitspreferentialattachmentinthetimecourseof RPI database accumulation. 

We selected RNA and protein complexes from the PDB database from – 2014to2023toconstructthedynamicprocess(seeSupplementaryFig.8a j) of the RNA–protein structural interaction network. Supplementary Fig. 9 shows the annual increase in nodes and edges within the network and displays the total number of nodes and edges. 

In the Methods section, we introduce a function that evaluates the attachment ability of nodes with varying degrees in the RNA–protein structural interaction network. This function measures the nodes’ potential to acquire new links over time. This capacity was computed for each period between 2014 and 2023. Afterward, we averaged the results across time intervals and refined them using Savitzky-Golay filtering (see Supplementary Fig. 10). 

Basedonthe findings,nodeswithdegreesrangingfrom38to55showa solid ability togainnewlinks.The abilityof RPInetworknodestoreach new nodes is directly related to their respective degrees. This result implies that RNA and protein nodes are preferred in the evolution of RPIs. Therefore, network features can be explored for deep learning in predicting RPIs. 

## Predicting RPI on benchmark datasets 

We estimated the prediction performance of ZHMolGraph on the benchmark datasets NPInter2 and RPI7317 using a 5-fold cross-validation. The 

Communications Biology | (2025) 8:247 

4 

Article 

https://doi.org/10.1038/s42003-025-07694-9 

dataset is partitioned into five segments, with four segments designated as the training set and the remaining one as the testing set. For each fold, the training set is further divided into ten parts. A random selection process designatesoneofthesepartsasavalidationsettofacilitatethedetermination of the optimal threshold for a binary decision of interactions. We determine the optimal threshold for the binary probability by selecting the threshold that maximizes the MCC value on the validation set. We compared ZHMolGraph with other best methods, including RPITER, RPISeq, IPMiner, and NPI-GNN. We calculated various measures, accuracy, sensitivity, specificity, precision, and MCC values, to evaluate the performance ofourmethod.Allevaluationmeasuresfordifferentmethodswereextracted from the ref. 23. 

As shown in Fig. 4a, ZHMolGraph achieved high accuracy, sensitivity, precision, specificity, and MCC with scoresof 0.955, 0.975, 0.935, 0.938,and 0.911, respectively. ZHMolGraph had the highest performance in MCC comparedtoother methods, making it a state-of-the-art method. To further validate the performance of ZHMolGraph, we conducted another 5-fold cross-validation on the RPI7317 dataset. ZHMolGraph achieved an accuracyof 0.915, a sensitivity of 0.927, a precision of 0.904,a specificity of 0.906, and an MCC of 0.831 (see Fig. 4b). ZHMolGraph still outperformed other methods with better MCC results. 

RNAs and proteins. ZHMolGraph achieved the best performance on unknown nodes, with average AUROC and AUPRC scores of 0.798 and 0.820, achieving improvements of 7.1% and 4.6% than RPITER in AUROC and AUPRC, respectively. Compared with IPMiner, we achieved13.4% and 8.1% enhancements. In addition, ZHMolGraph outperformed NPI-GNN with improvements of 28.7% and 30.0% in AUROC and AUPRC. We also tested the ability to identify correct interactions between ZHMolGraph and other methods. Figure 5b illustrates precision computed from the top 200, 600, 1000, 1400, 1800, and 2200 prediction probabilities. For ZHMolGraph, the average precision of the top 600 (1000) prediction probabilities stands at 0.854 (0.764), outperforming PRITER, IPMiner, and NPI-GNN, which report 0.826 (0.677), 0.814 (0.641) and 0.540 (0.529), respectively. The results suggest that most predictions of ZHMolGraph on unknown nodes are more reliable than other methods. 

The result of the TheNovel dataset proved that the deep learning method of graph neural networks heavily relies on known interaction network, making it challenging to apply to unknown nodes. The comparison between ZHMolGraph and models such as NPI-GNN and IPMiner confirms that unsupervised LLMs molecular embeddings enhance the generalizability of interaction prediction models. 

## The contribution of the LLMs and graph embeddings 

## Predicting RPI for unknown nodes 

To investigate the performance of ZHMolGraph further, we tested our method on even more challenging distinct scenarios, unknown nodes, when both RNAs and proteins from the test dataset are absent in the training data. We previously performed five-fold cross-validation training on the NPInter2 dataset, yielding five models. These models were then employed to test the unknown nodes testing set TheNovel, and the average performance of the five models was calculated to represent the final results. 

Figure 5ashowsthatthegraphneuralnetwork-basedmodelNPI-GNN performs poorly in the unknown testing, with average AUROCof 0.511 and AUPRC of 0.520. The deep learning-based method IPMiner shows performancewithaverageAUROCof0.664andAUPRCof0.739,andRPITER obtained average AUROC of 0.727 and AUPRC of 0.774 for unknown 

Large language models (LLMs) generate node embeddings that capture inter-residue dependencies to characterize molecules. In contrast, graph neural networks (GNNs) focus on the local environment within RNA–protein interaction networks, creating node embeddings incorporating information from known RNA–protein interactions. Combining these two types of information can enhance performance on existing net’ – works and improve the model s generalizability in predicting RNA protein interactions involving RNA and proteins that were not included during model training. 

We conducted a comparative analysis to evaluate the impact of large language models (LLMs) and graph neural network (GNN) embeddings. In our LLM approach, we utilized k-mer embedding methods that do not dependonLLMs,likethoseimplementedinNPI-GNNandIPMiner(please refer to the “Methods” section for more information). For the GNNs, we 

Fig. 4 | Overall performance of ZHMolGraph on the benchmark datasets. a, b The average performance comparison of ZHMolGraph, RPITER, NPIGNN, IPMiner, and RPISeq-RF on dataset NPInter2 (a) and PRI7317 (b) by 5-fold cross-validation. Source data are provided in Supplementary Data 6. 

Communications Biology | (2025) 8:247 

5 

Article 

https://doi.org/10.1038/s42003-025-07694-9 

Fig. 5 | Overall performance of ZHMolGraph on the unknown dataset TheNovel. a The average AUROC and AUPRC performance of ZHMolGraph, RPITER, IPMiner, and NPI-GNN on the unknown nodes dataset TheNovel, using n = 5 models trained on NPInter2 dataset (dots represent the performance of each model, bar height corresponds to the mean). b Top N prediction of ZHMolGraph, RPITER, IPMiner, and NPI-GNN on the unknown nodes dataset TheNovel, using n = 5 models trained on NPInter2 dataset (dots represent the precision of each model, the lines indicate the average precision for top-N predictions, and the shaded area around the average represents the standard deviation of the percision). Source data are provided in Supplementary Data 6. 

compared the performance of ZHMolGraph with a basic feedforward neural network that exclusively used the LLMs embeddings. 

The results from the benchmark dataset NPInter2 indicate that the LLMs method consistently outperforms methods without LLMs, achieving a higher average MCC value of 0.911 compared to 0.884 (see Fig. 6a). Similarly, in the benchmark dataset RPI7317, the LLMs method also attains a superior average MCC value of 0.831 compared to 0.815 (see Fig. 6b). Supplementary Figs. 11a and 11b show that including GNN embeddings improves the model’s accuracy by 2.4% and 3.6% on the NPInter2 and RPI7317 sets, respectively. The performance of the simplified model demonstrates that the graph structure indeed provides improvements. This suggests that the embeddings from LLMs and GNNs can complement each other, compensating for the respective limitations of each approach on benchmark tests. 

When it comes to the unknown testing dataset TheNovel, our method without LLMs performs poorly with 0.579 of average AUROC and 0.669 of average AUPRC. The results demonstrate that utilizing LLMs models of either protein or RNA leads to improvements. Specifically, the RNA LLMs model achieved an average AUROC of 0.738 and an average AUPRC of 0.770, while the protein LLMs model achieved an average AUROC of 0.636 and an average AUPRC of 0.707. Furthermore, the optimal performance is observed when employing both RNA and protein LLMs, resulting in an average AUROC of 0.798 and an average AUPRC of 0.820 (see Fig. 6c). We also evaluated the performance of ZHMolGraph by calculating the top 200, 600, 1000, 1400, 1800, and 2200 confident prediction probabilities for methods with and without LLMs models. Our method utilizing LLMs models achieves an average precision of 0.854 (0.764) for the top 600 (1000) confidentpredictionprobabilities.Incomparison,methodsrelyingsolelyon RNA-LLMs, protein-LLMs, and no-LLMs approaches yield precisions of 0.827(0.712),0.777(0.639),and 0.714(0.575),respectively(seeFig. 6d).The case-by-case comparison of prediction probabilities for positive samples between LLMs and no-LLMs methods reveals a superior performance by the LLMs methods (988/194), demonstrating the effectiveness of LLMs models (see Fig. 6e). The violin plots in Fig. 6f show the prediction probability of positive and negative samples from various methods with or without LLMs. Using LLMs methods enhances the distinguishability of prediction probabilities between positive and negative samples. The results (see Supplementary Fig. 11c) showed that the graph embedding method led to a 2% increase in AUROC and AUPRC regarding the unknown testing dataset TheNovel. 

These results indicate that when faced with entirely unknown RNA and protein sequences, the information from GNNs is limited to training the visible RPI network and provides limited assistance for unknown data. In contrast, LLMs exhibit stronger generalizability due to their advantages in understanding biomolecules, enabling ZHMolGraph to adapt to unknown nodes effectively. Detailed information about the 

contribution of the LLMs and graph embeddings can be found in Supplementary Tables 1 and 2. 

## Overcome the challenges of “orphan” RNAs and proteins 

“Orphan” nodes often lack direct connections or homologs to other nodes, which makes it very difficult to learn the features from such limited connections. We analyze the impact of “orphan” nodes in the benchmark dataset NPInter2 and RPI7317 concerning the prediction of interactions under the following scenarios: (i) unknown edges, where both RNA and protein from the testing set are present in training set folds during the fivefoldcross-validation;(ii)unknownRNA,whentheRNAfromthetestingset fold is absent in training set folds during the five-fold cross-validation; (iii) unknown protein, when the protein from the testing set fold is absent in training set folds during the five-fold cross-validation; and (iv) unknown nodes,whenboththeRNAandproteinfromthetestingsetfoldareabsentin training set folds during the five-fold cross-validation. In case (i), we perform random five-fold partitioning of the benchmark dataset. In cases (ii) and (iii), the RNA and protein nodes are equally divided into five groups, following the five-fold partitioning. The interactions corresponding to the RNA nodes in each group are designated as the testing set, while the remaining interactions serve as the training set. In case (iv), we utilize the training data from case (iii) and select interactions from the testing data where both the RNA and protein are not present in the training data for evaluation. As shown in Supplementary Fig. 12, the performance achieved using LLM embedding methods across the four scenarios (average MCC values of (i) 0.911/0.831, (ii) 0.859/0.778, (iii) 0.741/0.506, (iv) 0.549/0.515 for datasets NPInter2/RPI7317) is superior to that obtained without LLM embeddings (average MCC values of (i) 0.884/0.815, (ii) 0.825/0.744, (iii) 0.577/0.262, (iv) 0.362/0.129 for datasets NPInter2/RPI7317). 

An orphan degree (see Methods: The orphan degree definition for RNA/protein nodes from the testing set) for RNA/protein nodes from the testing set can be calculated. It can be observed that as the orphan degree increases, the methods utilizing LLMs demonstrate greater robustness compared to those without LLMs. In the case of an unknown protein, the RNA orphan degree is higher than the protein orphan degree in the unknown RNA scenario. This results in the model’s performance for unknown proteins being worse than that for unknown RNA (see Supplementary Fig. 12). These comparisons validate how LLMs molecular embeddings improve the generalizability of RPI prediction models. 

We further discuss the nodes with few or no homologs. It is generally accepted that an identity more significant than 35% indicates a similar structure for protein homologous sequences[29] , while the identity cutoff for RNA homologous sequences is 80%[30] . We compared the pairwise sequence similarity of nodes in NPInter2 (see Supplementary Figs. 13c, d for protein/ RNA) and created histograms to analyze the similarity distribution (see Supplementary Fig. 13b, e for protein/RNA). The results indicate that only a 

Communications Biology | (2025) 8:247 

6 

Article 

https://doi.org/10.1038/s42003-025-07694-9 

Fig. 6 | The contribution of the unsupervised LLMs. a, b The average performance comparison of ZHMolGraph on dataset NPInter2 and PRI7317 by 5-fold crossvalidation with or without LLMs methods (dots represent the performance of each fold, n = 5). c Our methods’ average AUROC and AUPRC performance on the unknown nodes dataset TheNovel, using n = 5 models trained with or without LLMs methods on the NPInter2 dataset (dots represent the performance of each model, bar height corresponds to the mean). d Top N prediction on the unknown dataset TheNovel, using n = 5 models trained with or without LLMs methods on NPInter2 dataset (lines indicate the average precision for top-N predictions, the shaded area 

tiny percentage (<1.6%) of node pairs are homologs in the two datasets. Besides, we analyze the effect of sequence similarity in the NPInter2 set. We calculated the average similarity of each protein node to all other nodes. Based on the sorted order of similarity, we performed five-fold partitioning. To ensure the discontinuity of similarity and the adequacy of training samples, we selected protein nodes corresponding to the top 0–10%, 20–30%, 40–50%, 60–70%, and 80–90% of sequence similarity as the test sets. At the same time, the remaining RPIs were used as the training set for five-fold cross-validation. As shown in Supplementary Fig. 13a, the LLM embeddings outperform the model without LLMs, demonstrating greater robustness. The Pearson correlation between MCC and similarity shows − 0.04 for LLMs and 0.10 for methods without LLMs (see Supplementary Fig. 13f). The methods without the LLMs embeddings exhibit a positive correlation with similarity. The decrease in performance for Fold 4 is attributed to a sudden increase in the orphan degree of RNA nodes in the Fold 4 test set (see Supplementary Fig. 13g). 

## Model performance across species 

We further conduct evaluate the bias of ZHMolGraph across different species. We trained the ZHMolGraph model on the benchmark NPInter2 and evaluated its performance across different species using TheNovel test set. The NPInter2 benchmark primarily includes three species: Homo sapiens, Mus musculus, and Saccharomyces cerevisiae. The proportions of 

around the average represents the standard deviation of the precision). e Comparison of the prediction probability for the positive samples between methods with or without LLMs. f Violin plots of the prediction probabilities for n = 1182 positive samples and n = 1182 negative samples, comparing methods with and without LLMs. In (a), (b) and (d), ZHMolGraph with both RNA and protein LLMs is colored in green, ZHMolGraph only with RNA LLMs is colored in yellow, ZHMolGraph only with protein LLMs is colored in purple, ZHMolGraph without LLMs is colored in red. Source data are provided in Supplementary Data 6. 

interaction edges for these species are 67.0%, 21.1%, and 8.6%, respectively, while the proportions of interaction nodes are 55.6%, 36.1%, and 4.2% (see Supplementary Fig. 14a, b). In TheNovel test set, the proportions of interaction edges for Homo sapiens and Mus musculus are 86.0% and 13.9%, respectively (see Supplementary Fig. 14c). 

As shown in Supplementary Fig. 14d, Mus musculus performs worse than Homosapiens inbothAUROCandAUPRC.Thisperformancefmight be due to the significantly lower number of Mus musculus samples in the training data NPInter2 compared to Homo sapiens. To investigate how data quantity affects model performance across different species, we conducted additional analysis by controlling the interaction ratios at 0.2, 0.4, 0.6, 0.8, and 1 for individual species, including Homo sapiens and Mus musculus, withintheNPInter2dataset.Weutilizedtrainingdataofvaryingsizesacross the species to train the models. Subsequently, we tested the models on TheNovel dataset. Supplementary Fig. 15a, b illustrates the model’s performanceundervaryingratiosof Homosapiens interactionsintheNPInter2 dataset, while Supplementary Fig. 15c, d presents the results for Mus musculus interactions. As the ratio of interactions for each species decreases, the model’s performance experiences a slight decline. However, the average AUROC and AUPRC remain around 0.7, highlighting the robustness of ZHMolGraph in maintaining accuracy across different species. When the ratio of Homo sapiens interactions is reduced to 20%, ZHMolGraph achieves an average AUROC of 0.705 and an AUPRC of 0.747 for Homo 

Communications Biology | (2025) 8:247 

7 

Article 

https://doi.org/10.1038/s42003-025-07694-9 

sapiens, and an average AUROC of 0.711 and an AUPRC of 0.772 for Mus musculus (see Supplementary Fig. 15a, b). Similarly, when the ratio of Mus musculus interactions drops to 20%, ZHMolGraph still delivers an average AUROC of 0.681 and an AUPRC of 0.743 for Homo sapiens, as well as an average AUROC of 0.661 and an AUPRC of 0.724 for Mus musculus (see Supplementary Fig. 15c, d). While a larger quantity of species can enhance the performance of corresponding species interactions, these results suggest that ZHMolGraph effectively minimizes potential structural biases toward specific RNA and protein families. 

We also investigate the effect of RPI network topological differences between species on model performance. As shown in Supplementary Fig. 16a, Homo sapiens has a higher proportion of hub nodes compared to Mus musculus, with an average degree of 4.901 (2.378 for Mus musculus). This indicates that each node in Homo sapiens carries more interaction information, contributing to its superior performance on the TheNovel test set compared to Mus musculus. The modularity of Homo sapiens and Mus musculus nodes is also investigated. As shown in Supplementary Fig. 16b, the Spearman correlation coefficients between the topological coefficients and degree distributions are −0.979 and −0.996 for Homo sapiens and Mus musculus,respectively.These resultsindicate thatboth species exhibit ahigh degree of modularity. The analysis of Homo sapiens and Mus musculus revealed that both species exhibit scale-free degree distribution and a highly modular structure. However, the differences in degree distribution between the two species might influence the model’s performance across the species. 

We further investigate the sequence similarity between species. Specifically, we examined the protein similarity between intra-species (see Supplementary Fig. 17a, b) and inter-species (see Supplementary Fig. 17c) for the two major species, Homo sapiens and Mus musculus. The histograms – (see Supplementary Fig. 17d f) reveal that only a few pairs of nodes (~1%) possess homologs (35% similarity) in intra-species and inter-species protein similarity comparisons, indicating that protein sequence similarity is not a primary influencing factor of ZHMolGraph performance. 

## Model performance with various data quality 

Since ZHMolGraph relies on the sequence data, we investigate the model’s performance if the input sequences contain errors or are incomplete. To introduce errors into the input sequences, we randomly masked residues/ nucleotides by replacing them with a ‘-’. The masking percentages are set at 1%, 2%, 3%, 4%, 5%, 10%, 15%, 20%, 25%, 30%, 35%, 40%, 45%, and 50%. We evaluated the performance of different mask percentages on the TheNovel test set. As shown in Supplementary Fig. 18, performance decreases by ~0.1 with a masking percentage of 1%. When the masking percentage is below 10%, the performance remains around 0.7 for the average AUROC and AUPRC. These results indicate that ZHMolGraph is robust, maintaining high accuracy with less than 10% sequence errors. Even when the sequence errors reach 50%, ZHMolGraph achievesan average AUROC and AUPRC of around 0.6. 

## Accurate identification of viral RNA–protein interactions 

We assess ZHMolGraph’s performance in identifying unknown viral RNA–protein interactions related to SARS-CoV-2 and compare it with IPMiner and RPITER. ZHMolGraph successfully identified 704 out of 819 viral interactions, achieving a recall rate of 0.860 (see Supplementary Fig. 19a). In contrast, IPMiner identified only 423 out of 819 viral interactions, resulting in a recall rate 0.518 (see Supplementary Fig. 19b). In comparison, RPITER identified 68 out of 819viral interactions, with a recall rate of 0.083 (see Supplementary Fig. 19c). These results demonstrate ZHMolGraph’s superior potential for uncovering viral RNA–protein interactions. 

## Discussion 

The precise identity of RNA–protein interactions is a prerequisite for RNA–protein complex prediction. The large language models interpret an entire sequence as one sentence and its amino acids (nucleotides) as individual words. We designed ZHMolGraph to improve accuracy and 

generalizability by combining unsupervised LLMs language models and graphneuralnetworkstoenhanceaccuracy.Asa result,ourmethodshowed a noticeable performance improvement compared to other deep learning methods, with a commendable increase in AUROC (13.4%) and AUPRC (8.1%). 

– We can utilize ZHMolGraph to help predict the RNA protein complex structure. ZHMolGraph can identify the sequences in the interface where RNAand protein bind together. To locate the binding spots along the sequence, we perturb the nucleotides and observe the changes in ZHMolGraph prediction. We select the tenlowest resulting binding probabilities as the binding locations along the sequence (see Methods: Utilizing ZHMolGraph for RNA-protein structure prediction). 

Then, we applied the location information to address a challenging scenario of fully flexible unbound RNA–protein complex prediction. We tested on two unbound-unbound RNA–protein complexes and two unbound-bound RNA–protein complexes. We examined the predicted – RNAbinding site locations in the four RNA protein complexes. The results have shown that an average of 70% of the predicted RNA binding sites corresponded to the correct ones. In contrast, the remaining incorrectly predicted binding sites were mainly located near the native binding sites, demonstrating the accuracy of the interface recognition (see Supplementary Fig. 20). 

For the 1JID structure, the lowest RMSD model from 3dRPCZHMolGraph is 2.64 Å, compared to 4.82 Å for 3dRPC, with a TM-score of 0.83 versus 0.79. 3dRPC-ZHMolGraph achieved an average RMSD of 9.60 Å (colored in red), outperforming RMSDs of around 11.01 Å in the 3dRPC (colored in blue) (see Fig. 7a). In the structure of 3HL2, the lowest RMSD model of 3dRPC-ZHMolGraph is 3.71 Å, compared with 5.86 Å for 3dRPC, with a TM-score of 0.91 versus 0.88. 3dRPC-ZHMolGraph achieved an average RMSD of 13.19 Å (colored in red), better than RMSDs of around 19.21 Å in the 3dRPC (colored in blue) (see Fig. 7b). In the structure of 2ANR, the lowest RMSD model of 3dRPC-ZHMolGraph is 0.94 Å, compared with 5.83 Å for 3dRPC, with a TM-score of 0.98 versus 0.88. 3dRPC-ZHMolGraph achieved an average RMSD of 9.07 Å (colored in red), better than RMSDs of around 9.46 Å in the 3dRPC (colored in blue) (see Fig. 7c). In the structure 1N78, the lowest RMSD model of 3dRPC-ZHMolGraph is 1.92 Å, compared with 9.73 Å for 3dRPC, with a TM-score of 0.94 versus 0.84. 3dRPC-ZHMolGraph achieved an average RMSD of 15.32 Å (colored in red), better than RMSDs of around 18.78 Å in the 3dRPC (colored in blue) (see Fig. 7d). Together, 3dRPC-ZHMolGraph led to a shift in the distribution of the predictions toward lower RMSDs than 3dRPC (Histograms in Fig. 7). These findings suggest that the binding information provided by ZHMolGraph can enhance the performance of RNA–protein structure – prediction. The detail of the unbound RNA protein complex prediction performance can be found in Supplementary Table 3. 

Although ZHMolGraph can contribute to identifying RPI interactions with sequence information, the binding phenomenonprimarily depends on the tertiary structures of molecules. These structuresshape the architectures of binding pockets and influence bond rotations[31] . As the number of experimentally determined 3D structures of RNA–protein complexes continues to increase and be integrated, it will significantly enhance our research on ZHMolGraph. By incorporating higher-order molecular properties that drive RNA–protein binding, this progress will further improve the predictive accuracy of ZHMolGraph, allowing for a more detailed analysis of physical interaction contacts between RNA–protein residues and nucleotides. 

In conclusion, we have comprehensively analyzed distinct features within the RPI interaction network from diverse perspectives and proposed an innovative AI-driven approach to predicting RNA–protein interactions, integrating large language models and graph neural networks. Benchmark evaluations on two benchmark datasets demonstrate that ZHMolGraph exhibits high accuracy and generalizability compared to other state-of-the-art methods, even for unknown node datasets. ZHMolGraph is a reliable tool for determining RNA–protein 

Communications Biology | (2025) 8:247 

8 

Article 

https://doi.org/10.1038/s42003-025-07694-9 

– Fig. 7 | The unbound RNA protein complex prediction ranking analysis. a Human SRP19 in complex with helix 6 of human SRP RNA (PDB code: 1JID). b Human SepSecS-tRNASec complex (PDB code: 3HL2). c Nova-1 KH1/KH2 domain tandem with RNA hairpin (PDB code: 2ANR). d Glutamyl-tRNA synthetase complexed with tRNA (Glu) and glutamol-AMP (PDB code: 1N78). Histograms (left to right) are the RMSD (Root Mean Square Deviation) distributions of the 3dRPC-ZHMolGraph with binding information (red) compared to the 3dRPC 

interactions on a genome-wide scale and accurately predicting nearnative RNA–protein structures. We expect the method to be helpful for RNA–protein complex prediction and drug design. 

## Methods 

## Statistics on RPI network data from various repositories 

To enhance our comprehension of the biological processes controlled by RPIs, we utilized three distinct sources of experimental data: structural data from the Protein Data Bank (PDB), high-throughput data from the RNAInter database, and literature mining validated data from the NPInter5 database[32][–][34] . Our approach considers RNA/proteins as individual nodes and their interactions as edges, enabling us to construct a molecular-level network. 

without binding information (blue), relative to the native structure among the top 50 predictions. When binding information is considered, the 3dRPC-ZHMolGraph distribution shifts toward lower RMSD values. The lowest RMSD models in the top 50 predictions of 3dRPC-ZHMolGraph (RNA in red) are more similar to native complexes (RNA in green) than 3dRPC (RNA in blue). Source data are provided in Supplementary Data 6. 

## Structural interactions 

We extracted the RNA–protein structural interactions from RNA–protein complex structures in the Protein Data Bank (before October 28, 2023) throughthefollowingprocess[32] :(i)Westartedbydownloadingallstructures with an entry polymer type labeled ‘Protein/NA’. We kept the first model from each of these structures. (ii) We then filtered out all DNA structures. (iii) We removed RNA chains containing non-AUCG bases in their sequences. (iv) RNA and protein chains with less than 20 units were excluded. As a result, we compiled a set of PDB structures that exclusively contained RNA–protein complexes. Using these structures, we constructed structure-basedRPInetworks.Intheresidue/nucleotide-levelnetwork,each node symbolizes a residue/nucleotide, and an edge signifies an interaction between the residue and nucleotide when the non-covalent distance 

Communications Biology | (2025) 8:247 

9 

Article 

https://doi.org/10.1038/s42003-025-07694-9 

between any pairof their heavy atomsisless than 8 Å. Inthe molecular-level network, nodes represent RNA and protein, and edges represent interactions between RNA and protein if an interaction occurs between residues/ nucleotides within the RNA and protein nodes. In addition, RNA and protein chains having the same sequence are treated as identical nodes. 

## High-throughput interactions 

RNAInter is a database consolidating RNA interaction data validated by experiments and obtained from scientific literature and databases. The interactions are matched with comprehensive RNA annotations, including RNA editing, localization, modifications, and associated diseases[33] . Moreover, it provides a confidence scoringsystemthat assessesthe quality of each interaction. This is done by integrating the reliability of experimental evidence, scientific consensus, and types of tissues/cells involved. We specifically extracted interactions involving Homo sapiens from the RNAInter database. These interactions were determined using the high-throughput technique HITS-CLIP and then filtered according to the recommended confidencescorecutoff of 0.186,asmentioned intheRNAInterdatabase.As a result, we obtained a dataset comprising 82,170 high-throughput interactions involving 23,205 RNAs and 4713 proteins, which was used to construct the high-throughput interaction network. 

## Literature-mining validated interactions 

The NPInter5 database comprehensively collects interactions between noncoding RNAs and various biomolecules such as proteins, RNAs, and DNAs[34] . The interactions are from different sources, such as experimentally validateddatafromliteratureminingandhigh-throughputdatasets.Forthis study,wefocused onRPIssupported byexperimentalevidencedocumented in published literature. We extracted the interactions from the NPInter database to form the literature-mining validated RPI network, resulting in a collection of 4329 interactions involving 1937 RNAs and 1716 proteins. 

## Topological coefficient calculation 

The topological coefficient measures the degree to which a node shares its neighbors with othernodesinthenetwork.The topologicalcoefficientvalue TC of the node v can be calculated using the following formula: 

**==> picture [166 x 22] intentionally omitted <==**

Inthisequation, m representsa nodethat shares atleastoneimmediate neighbor with v. J(v, m) denotes the count of shared neighbors between v and m, plus one if v and m are connected directly[35] . 

## Attachment rate calculation 

The attachment rate measures nodes’ varying ability to acquire new connections. Forinstance,suppose we haveinformationregardingthestructure of a network at two different time points. The first time point is marked as T0, while the subsequent time point is denoted as T1, which is T1 = T0 + ΔT. We use the function Γ(k) to represent the capacity of a node with degree k to obtain new connections, which is determined using the following formula: 

**==> picture [223 x 26] intentionally omitted <==**

where Nk stands for the count of nodes with a degree of k at the time T0, and degree(n,T0) representsthedegreeofnode n atthetime T0.Thefunction Γ(k) characterizes the ability of nodes with degree k from T0 to T1 to acquire new attachments. A stronger (weaker) Γ(k) value indicates that nodes with degree k[36] possess a more significant (lesser) capacity to establish new connections during the dynamic process. 

dimensional embeddings for RNAs and 1024-dimensional embeddings for proteins. The LLMs models for the embeddings are unsupervised, meaning they are conducted independently of the subsequent interaction prediction task. Given a nucleotide sequence x0 and an amino acid sequence x1, we encode them using RNA-FM and ProtTrans: 

**==> picture [226 x 14] intentionally omitted <==**

**==> picture [171 x 10] intentionally omitted <==**

Graph neural network embedding layer. We use an inductive node representation learning graph neural network, which offers an efficient approach to generating node representations from the available RPI network[24] . The embeddings from LLMs models are passed into the graph neural network module. We pad the RNA embedding vector with zeros to make it a 1024-dimensional vector. We used two SAGE layers as our graph neural network module, which can rapidly learn a function that produces local embeddings by sampling and aggregating features from nodes’ local neighborhoods without requiring knowledge of the Laplacian matrix of the entire network. The module randomly samples neighbors from a target node’s first-order neighborhood. Then, for each sampled neighbor, samples from their neighborhood are recursively taken in the following layer. This process reduces computational complexity while preserving multi-hop neighbor information. After each sampling layer, the module aggregates the feature vectors of neighboring nodes using a learnable aggregation function. The aggregation considers the features of the neighboring nodes and their relative relationships, ensuring that the generated node embeddings have strong generalization capabilities. Ultimately, we obtain a fixed-dimensional embedding representation for each node in the RNA–protein interaction graph. More formally, for each node v with LLMs embedding �xðv[0] Þ in the RPI network, two SAGE layers compute: 

**==> picture [215 x 45] intentionally omitted <==**

**==> picture [153 x 30] intentionally omitted <==**

where k = 1, 2. N(v) denotes the neighborhood of node v,||∙|| represents the Euclidean norm. zv ¼ �xv[ð][2][Þ][denotes the node][ v][ output representation of the] graph neural network module. Various aggregator architectures can aggregate the neighbor representations, including the mean and max poolingaggregators.Acomparisonof differentaggregationfunctionscanbe found in Supplementary Note 1 and Supplementary Fig. 21. To make the aggregator function trainable and maintain a high representation of the neighborhood information, we used the mean aggregator in ZHMolGraph in our work. 

The graph neural network leverages node features of node profile information and node degrees to learn an embedding function that generalizes to unknown nodes. When testing the unknown nodes, we use the trained graph neural network model to generate embeddings by applying the learned aggregation functions. We set the out_channels of the graph neural network to be 100 and employed the unsupervised loss as the loss function. This loss function encourages similar representations for nearby nodes while ensuring that distant nodes remain distinct. For each output representations zu, 8u ∈ V, the graph-based loss function is defined as: 

## ZHMolGraph deep learning architecture 

LLMs embedding layer. ZHMolGraph utilizes LLMs RNA-FM and ProtTrans models (see Fig. 2b)[27][,][28] . These models generate 640- 

**==> picture [232 x 18] intentionally omitted <==**

Communications Biology | (2025) 8:247 

10 

Article 

https://doi.org/10.1038/s42003-025-07694-9 

where v is a node that co-occurs near u, σ is the sigmoid function, Pn is a negative sampling distribution, and Q defines the number of negative samples, which is set to 10. The first part of the loss function is designed to ensure that if nodes u and v are close together in the actual graph, their node embeddings should be semantically similar. Conversely, the second part of the loss function aims to achieve the opposite: if nodes u and v are distant in the actual graph, we expect their node embeddings to differ. In addition, we present the convergence curve of the graph neural network’s loss after 10 epochs (see Supplementary Figs. 22 and 23), which illustrates that the training process has converged successfully. When the graph neural network training was completed, the embedding vectors were subsequently fed into the VecNN training phase. 

VecNN architecture. We concatenate LLMs embedding features with graph neural network embedding features and then feed them into a neural network architecture, VecNN. The calculations follow these formulas: 

**==> picture [214 x 35] intentionally omitted <==**

**==> picture [204 x 108] intentionally omitted <==**

We chose a 1D convolutional layer as the feature extraction layer (a detailed comparison between different feature extraction layers can be found in Supplementary Note 2 and Supplementary Fig. 24), and the concatenated RNA and protein embeddings are individually input into the 1D convolutional layer with a kernel size of 3 to extract their features. These features then go through a dense layer that outputs a 2048-dimensional vector to standardize the dimensions. Next, the RNA and protein embeddings are combined and sent through two fully connected dense layers, resulting in the final binding probability output. 

## Model training and testing 

We trained VecNN by minimizing the loss function of binary cross-entropy usingtheAdamoptimizerwitha learningrate lr =0.001, β1 =0.9, β2 =0.999, abatchsizeof128,andanepochnumberof120.Acosineannealinglearning rate scheduler was employed. The hyperparameters optimization and the loss convergence curve can be found in Supplementary Note 3 and Supplementary Figs. 25–27. An analysis of the model’s computational complexity and scalability regarding the time and GPU memory requirements are provided in Supplementary Note 4 and Supplementary Fig. 28. 

Weemploythesamerandom five-foldcross-validationstrategyusedin previous studies to train the VecNN in the benchmark sets to ensure accuracy. First, a random five-fold split on the benchmark datasets can ensure balanced distributions of RNA–protein pairs in each fold. First, a random five-fold split ensures a balanced distribution between positive and negative samples. The ratios of positive to negative samples in the training sets were 1.000, 0.996, 0.992, 0.999, and 0.997, respectively. Second, the random five-fold partitioning has the lowest orphan degree and the best performance compared to other partitions, demonstrating that random five-fold partitioning possesses a balanced distribution of nodes (see Results: Overcome the challenges of “orphan” RNAs and proteins). 

We employ the five-fold cross-validation models trained on the NPInter2 dataset to evaluate the performance of ZHMolGraph on the 

dataset TheNovel and the viral RPI related to SARS-CoV-2. We set the output probability threshold to 0.5 for ZHMolGraph and other comparable methods to calculate the recall of viral RPI identification. 

## The large language models for nodes embedding 

With the rapid advancement of sequencing technology, there has been an exponential accumulation of unannotated RNA and protein sequences. Conversely, experimental annotation and the determination of 3D structures are expensive and yield only limited datasets. Consequently, we have utilized two unsupervised large language models, RNA-FM and ProtTrans, to encode genome-based features of RNA and proteins, respectively[27][,][28] . 

RNA-FM. The RNA-FM model is constructed upon 12 transformer-based bidirectional encoderblocks and is trained using self-supervised learning on 23 million ncRNA sequences sourced from the RNAcentral database. This self-supervised learning approach is based on the architecture of the BERT language model. RNA-FM acquires an understanding of sequential distribution and patternsthat encapsulate underlying structural and functional information through this process. It has been demonstrated that embeddings generated by RNA-FM consistently outperform state-of-the-art approaches across various downstreamprediction tasks related to structural and functional aspects, including RNA secondary structure prediction and 3D modeling results. For each RNA sequence of length L inputted into RNA-FM, it generates an embedding matrix of size L × 640. We subsequently compute the average along the L dimension of the L × 640 embedding matrix to derive a one-dimensional embedding vector of size 640. Given RNA-FM’s 1024-nucleotide input limit per RNA sequence, we process sequences that exceed this length as follows: (i) we divide them into segments of 1024 nucleotides; (ii) each segment is then processed sequentially by RNA-FM; and (iii) we average the embeddings from each segment to generate a single, one-dimensional embedding vector of size 640. This approach enables us to process longer sequences while retaining the overall context. Our previous research has demonstrated that the binding fragments in RNA–protein structures typically span approximately 10 nucleotides, much shorter than the 1024-nucleotide limit[37] . Therefore, the impact of segmenting sequences on critical motifs or interactions that span across chunk boundaries is likely negligible. 

ProtTrans. The ProtTrans models are language models designed for proteins that train two auto-regressive models (Transformer-XL and XLNet) and four auto-encoder models (BERT, Albert, Electra, and T5). These models were trained on data from UniRef and BFD, which collectively encompass up to 393 billion amino acids. Among these language models, the most informative embeddings model ProtT5-XL-UniRef50, utilizing an auto-encoder T5 architecture trained with 8-way model parallelism and 3 billion parameters on the UniRef50 dataset, has shown excellent performanceindownstreamtaskssuchasproteinsecondarystructurepredictions. It has even outperformedESM-1b,anunsupervised protein languagemodel based on Transformer architecture, which has been successfully employed in protein structure prediction. We choose the ProtTrans model ProtT5XL-UniRef50 as our protein LLMs embedding model. For every protein sequence with a length of L inputted into ProtTrans, the model produces an embedding matrix sized L × 1024. Following this, we calculate the average along the L dimension of the L × 1024 embedding matrix to obtain a onedimensional embedding vector sized 1024. Since there are other LLMs like ESM-2[38] forproteins,acomparativeanalysisbetweenProtTransand ESM-2 has been conducted and is detailed in Supplementary Note 5, Supplementary Fig. 29, and Supplementary Tables 4 and 5. 

## Datasets for model testing 

We collected two published benchmark datasets to test our model. These datasets are NPInter2 and RPI7317 datasets[36][,][39] . The NPInter2 dataset encompasses 10,412 experimentally confirmed NPIs, involving 4636 RNAs and449proteinsfromsixmodelspecies.ThesequencesoftheRNAsandthe proteins were sourced from the NONCODE database and the UniProt 

Communications Biology | (2025) 8:247 

11 

Article 

https://doi.org/10.1038/s42003-025-07694-9 

database, respectively[40][,][41] . On the other hand, the RPI7317 dataset contains 7317 experimentally validated NPIs involving 1874 RNAs and 118proteins. It was obtained from the NPInter3 database by selecting only interactions from humans. The sequences of RPI7317 were obtained from the Gencode and the UniProt databases[41][,][42] . As these two datasets exclusively contain positive samples, we generated negative samples of equal size by randomly selecting RNAs and proteins lacking validated interactions, following the same strategy used in IPMiner and NPI-GNN. To evaluate the impact of random negative sample selection on model performance, we conducted 10 trials using different random selections of negative samples and compared the results. As illustrated in Supplementary Fig. 30, the performance across various metrics is highly consistent. For the NPInter2 test set, the standard deviations for accuracy, recall, specificity, precision, and MCC are 0.002, 0.003, 0.002, 0.002, and 0.003, respectively (see Supplementary Fig. 30a). In the RPI7317 test set, the standard deviations for these five metrics are 0.004, 0.014, 0.008, 0.006, and 0.009, respectively (see Supplementary Fig. 30b). These findings suggest that the model exhibits strong robustness across different selections of negative samples. 

To evaluate the performance of ZHMolGraph on totally unknown RNAs and proteins, we further constructed an unknown testing set “TheNovel” from the NPInter5 dataset, which contains biologically significant RNA–protein experimental validated interactions obtained through literature mining. We excluded RNAs and proteins in the NPInter5 dataset, which were nodes in the NPInter2 dataset. This precaution ensured that the TheNovel testing set did not overlap with the NPInter2 dataset. After that, we downloaded the sequences of RNA and protein from the NONCODE database and the UniProt database, respectively[40][,][41] . We obtained 1182 interactions for unknown testing in the TheNovel dataset. 

To evaluate the ability of ZHMolGraph on identifying the viral RNA–protein interactions, we collected 1004 records of RNA–protein interactions related to SARS-CoV-2 from the NPInter5 database. After – removing duplicates, we had 819 unique RNA protein interactions. 

## Performance metrics 

We utilized accuracy, sensitivity, specificity, precision, and Mathews correlation coefficient (MCC) to evaluate the performance of model ZHMolGraph. We calculated the evaluation metrics as follows: 

**==> picture [188 x 20] intentionally omitted <==**

**==> picture [168 x 20] intentionally omitted <==**

**==> picture [168 x 20] intentionally omitted <==**

**==> picture [165 x 20] intentionally omitted <==**

**==> picture [235 x 24] intentionally omitted <==**

The TP (FP) represents true (false) positive, i.e., the number of true (false)interactions among the predictedinteractions. TheTN represents the true negative, i.e., the number of true non-interactions among the predicted non-interactions. The FN represents false negative, i.e., the number of true interactions the method fails to predict. Accuracy is the proportion of the correctlypredictedsamples.Sensitivity givesthe fractionof true interactions that can be predicted. Specificity gives the fraction of true non-interactions that can be predicted. Precision is the proportion of the correctly predicted interactions out of the total number of predicted interactions. MCC is a correlation coefficient between the true and predicted classes. 

## No-LLMs embedding methods 

In contrast to LLMs embedding methods, we utilized no-LLMs embedding methods for RNA and protein. RNA and protein sequences are encoded using the feature representation methods previously employed by Shen et al.[23] . For RNA, we utilized 3-mers, which generated 64dimensional embeddings for each RNA sequence. For protein, the 20 amino acids are classified into 7 groups based on their dipole moments and side chains volume: {A, G, V}, {I, L, F, P}, {Y, M, T, S}, {H, N, Q, W}, {R, K}, {D, E}, {C}. We employed 2-mers to encode protein embeddings according to these 7 groups, resulting in 49-dimensional embeddings for each protein sequence. 

## – Utilizing ZHMolGraph for RNA protein structure prediction 

ZHMolGraph can be used to find potential binding interfaces on the sequence. We plan to utilize this interface information to improve the RNA–protein structure prediction. Initially, We perturb each nucleotide in the RNA sequence by substituting it with ‘X’, then use ZHMolGraph to – calculatetheprobabilityofRNA proteininteractionsaftertheperturbation. The ten lowest resulting binding probabilities indicate the sequence’s most predictive RNA binding sites. 

Wethenuse3dRPC[15][,][43][,][44] ,atoolforRNA–proteinmoleculedocking,to generate 1000 RNA–protein unbound docking structures for these four examples. We then filter the 3dRPC docking results using binding information provided by ZHMolGraph. Specifically, within the 1000 RNA–protein docking structures, we retained a structure if the residues at the interface contained over 50% of the predicted top 10 active binding sites; otherwise, we discarded it. 

We analyzed the top 50 best-scoring structural models from 3dRPC and 3dRPC with ZHMolGraph-predicted binding information (3dRPC– ZHMolGraph). We tested on two unbound RNA protein complexes: the human SRP19 in complex with helix 6 of human SRP RNA (PDB code: 1JID)[45] and the human SepSecS-tRNASec complex (PDB code: 3HL2)[46] . In addition, we included two unbound RNA–protein complexes: Nova-1 KH1/KH2 domain tandem with RNA hairpin (PDB code: 2ANR) with unbound RNA and bound protein[47] , and Glutamyl-tRNA synthetase complexed with tRNA (Glu) and glutamol-AMP (PDB code: 1N78) with bound RNA and unbound protein[48] . 

For each RNA–protein complex, 1000 predictive structures were generatedby3dRPC,usingthecommandof ‘3dRPC-mode9-system8-par ’ RPDOCK.par . The CAPRI criterion evaluates the quality of the RNA–protein complexprediction. The Irmsd is the interface RMSD between the predicted and native structures after the superposition of corresponding proteins. In this study, ‘RMSD’ refers to the interface RMSD. RMSD calculations were performed by 3dRPC using the command of ‘3dRPC –mode 2 –system 0 –par RMSD.par’. In addition, we have incorporated the TMscore to evaluate the complex structures, which are calculated using the USalign tool[49] . 

## The orphan degree definition for RNA/protein nodes from the testing set 

We define an orphan degree for RNA/protein nodes from the testing set as follows: 

**==> picture [241 x 27] intentionally omitted <==**

where f ðxÞ ¼ 1þ1e[�][x][,][N][r][=][p][represents][the][number][of][RNA/protein][nodes] from the testing set, Gt represents the testing set RPI network, kr=p represents the degree of testing RNA/protein nodes in the training set network. The orphan degree represents the multiplicative inverse of the average degree of testing RNA/protein nodes within the training set network. A higher orphan degree indicates that the test nodes appear less frequently in the training set, further reflecting the degree of “ ” orphan . 

Communications Biology | (2025) 8:247 

12 

Article 

https://doi.org/10.1038/s42003-025-07694-9 

## Statistics and reproducibility 

Correlation analysis was conducted by calculating the Spearman and Pearson correlation coefficients using Python (version 3.8) with the scipy package(version1.10.1).Thepowerlawwithdegreeexponentsinthedegree distribution of RNA–protein interaction networks was performed using Python (version 3.8) with the powerlaw package (version 1.5). 

To ensure reproducibility, all information regarding model training and testing design can be found in the Methods and Supplementary Notes. The data supporting the experimental results can be found in the “Data availability” section, and the code supporting the results can be found in the “Code availability” section. The source data behind the graphs in the paper are listed in Supplementary Data 5 and 6. 

## Reporting summary 

Further information on research design is available in the Nature Portfolio Reporting Summary linked to this article. 

## Data availability 

The data that support the results of thisstudy are available from the authors. The data on the RNA–protein interaction networks are available in the Supplementary files. The datasets used to evaluate the ZHMolGraph are available at https://zenodo.org/records/14747845. Source data for Supplementary Fig.13d, ecan be found in SupplementaryData 5. Theother source data behind graphs in the paper are provide in Supplementary Data 6. 

## Code availability 

The code for validating the results of ZHMolGraph is provided at https:// github.com/Zhaolab-GitHub/ZHMolGraph[50] . 

## Received: 25 August 2024; Accepted: 6 February 2025; 

## References 

1. Cramer, P. Organization and regulation of gene transcription. Nature 573, 45–54 (2019). 

2. Zhao, B. S., Roundtree, I. A. & He, C. Post-transcriptional gene regulation by mRNA modifications. Nat. Rev. Mol. Cell Biol. 18, 31–42 (2017). 

3. Statello, L., Guo, C.-J., Chen, L.-L. & Huarte, M. Gene regulation by long non-coding RNAs and its biological functions. Nat. Rev. Mol. Cell Biol. 22, 96–118 (2021). 

4. Gebauer, F., Schwarzl, T., Valcárcel, J. & Hentze, M. W. RNA-binding proteins in human genetic disease. Nat. Rev. Genet. 22, 185–198 (2021). 

5. McCune, J. M. The dynamics of CD4+ T-cell depletion in HIV disease. Nature 410, 974–979 (2001). 

6. Wang, H. et al. A computational study of Tat-CDK9-Cyclin binding dynamics and its implication in transcription-dependent HIV latency. Phys. Chem. Chem. Phys. 22, 25474–25482 (2020). 

7. Szpotkowski, K., Wójcik, K. & Kurzyńska-Kokorniak, A. Structural studies – 

of protein nucleic acid complexes: a brief overview of the selected techniques. Comput. Struct. Biotechnol. J. 21, 2858–2872 (2023). 

8. Guo, L. et al. Biochemical and structural insights into RNA binding by ssh10b, a member of the highly conserved Sac10b protein family in archaea. J. Biol. Chem. 289, 1478–1490 (2013). 

9. Ahmed, M., Marchanka, A. & Carlomagno, T. Structure of a proteinRNA complex by solid-state NMR spectroscopy. Angew. Chem. 59, 6866–6873 (2020). 

10. Nakagawa, R. et al. Structure and engineering of the minimal type VI CRISPR-Cas13bt3. Mol. Cell 82, 3178–3192.e3175 (2022). 

11. Danan, C., Manickavel, S. & Hafner, M. PAR-CLIP: a method for transcriptome-wide identification of RNA binding protein interaction sites. Methods Mol. Biol. 1358, 153–173 (2016). 

12. Ray, D. et al. Rapid and systematic analysis of the RNA recognition specificities of RNA-binding proteins. Nat. Biotechnol. 27, 667–670 (2009). 

13. Jain, R. et al. RIP-Chip analysis: RNA-binding protein immunoprecipitation-microarray (Chip) profiling. Methods Mol. Biol. 703, 247–263 (2011). 

14. Darnell, R. B. HITS-CLIP: panoramic views of protein-RNA regulation in living cells. Wiley Interdiscip. Rev. RNA 1, 266–286 (2010). 

15. Huang, Y., Liu, S., Guo, D., Li, L. & Xiao, Y. A novel protocol for threedimensional structurepredictionof RNA-protein complexes. Sci.Rep. 3, 1887 (2013). 

16. Zheng, J., Hong, X., Xie, J., Tong, X. & Liu, S. P3DOCK: a protein-RNA docking webserver based on template-based and template-free docking. Bioinformatics 36, 96–103 (2020). 

17. Tuszynska, I., Magnus, M., Jonak, K., Dawson, W. & Bujnicki, J. M. NPDock: a web server for protein-nucleic acid docking. Nucleic Acids Res. 43, W425–W430 (2015). 

18. Liu, H., Zhuo, C., Gao, J., Zeng, C. & Zhao, Y. AI-integrated networkfor RNA complex structure and dynamic prediction. Biophys. Rev. 5, 041304 (2024). 

19. Berman, H. M. et al. The Protein Data Bank. Nucleic Acids Res. 28, 235–242 (2000). 

20. Abramson, J. et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. Nature 630, 493–500 (2024). 

21. Muppirala, U. K., Honavar, V. G. & Dobbs, D. Predicting RNA-protein interactions using only sequence information. BMC Bioinforma. 12, 489 (2011). 

22. Pan, X., Fan, Y.-X., Yan, J. & Shen, H.-B. IPMiner: hidden ncRNAprotein interaction sequential pattern mining with stacked autoencoder for accurate computational prediction. BMC Genomics 17, 582 (2016). 

23. Shen, Z. A., Luo, T., Zhou, Y. K., Yu, H. & Du, P. F. NPI-GNN:predicting ncRNA-protein interactions with deep graph neural networks. Brief. Bioinforma. 22, bbab051 (2021). 

24. Hamilton, W. L., Ying, R. & Leskovec, J. Inductive representation learning on large graphs. In: Proceedings of the 31st International Conference on Neural Information Processing Systems) (Curran Associates Inc., 2017). 

25. Chowdhury, R. et al. Single-sequence protein structure prediction using a language model and deep learning. Nat. Biotechnol. 40, 1617–1623 (2022). 

26. Singh, A. Structure prediction for orphan proteins. Nat. Methods 20, 176–176 (2023). 

27. Elnaggar, A. et al. ProtTrans: toward understanding the language of life through self-supervised learning. IEEE Trans. Pattern Anal. Mach. Intell. 44, 7112–7127 (2022). 

28. Chen, J. et al. Interpretable RNA foundation model from unannotated data for highly accurate RNA structure and function predictions. Preprint at bioRxiv https://doi.org/10.1101/2022.08.06.503062 (2022). 

29. Olivella, M., Gonzalez, A., Pardo, L. & Deupi, X. Relation between sequence and structure in membrane proteins. Bioinformatics 29, 1589–1592 (2013). 

30. Zhang, C., Zhang, Y. & Pyle, A. M. rMSA: a sequence search and alignment algorithm to improve RNA structure modeling. J. Mol. Biol. 435, 167904 (2023). 

31. Liu, H., Jian, Y., Hou, J., Zeng, C. & Zhao, Y. RNet: a network strategy to predict RNA binding preferences. Brief. Bioinform. 25, bbad482 (2023). 

32. Burley, S. K. et al. RCSB Protein Data Bank: biological macromolecular structures enabling research and education in fundamental biology, biomedicine, biotechnology and energy. Nucleic Acids Res. 47, D464–D474 (2019). 

33. Kang, J. et al. RNAInter v4.0: RNA interactome repository with redefined confidence scoring system and improved accessibility. Nucleic Acids Res. 50, D326–d332 (2022). 

34. Zheng, Y. et al. NPInter v5.0: ncRNAinteraction database in a new era. Nucleic Acids Res. 51, D232–d239 (2023). 

Communications Biology | (2025) 8:247 

13 

Article 

https://doi.org/10.1038/s42003-025-07694-9 

35. Stelzl, U. et al. A human protein-protein interaction network: a resource for annotating the proteome. Cell 122, 957–968 (2005). 

36. Fan, X.-N. & Zhang, S.-W. LPI-BLS: predicting lncRNA–protein interactions with a broad learning system-based stacked ensemble classifier. Neurocomputing 370, 88–93 (2019). 

37. Yang, R. et al. RPpocket: an RNA–protein intuitive database with RNA pocket topology resources. Int. J. Mol. Sci. 23, 6903 (2022). 

38. Lin, Z. et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science 379, 1123–1130 (2023). 

39. Yuan, J. et al. NPInter v2.0: an updated database of ncRNA interactions. Nucleic Acids Res. 42, D104–D108 (2014). 

40. Bu, D. et al. NONCODE v3.0: integrative annotation of long noncoding RNAs. Nucleic Acids Res. 40, D210–D215 (2012). 

41. Wang, Y. et al. A crowdsourcing open platform for literature curation in UniProt. PLOS Biol. 19, e3001464 (2021). 

42. Frankish, A. et al. GENCODE: reference annotation for the human and mouse genomes in 2023. Nucleic Acids Res. 51, D942–d949 (2023). 

43. Huang, Y., Li, H. & Xiao, Y. 3dRPC: a web server for 3D RNA–protein structure prediction. Bioinformatics 34, 1238–1240 (2017). 

44. Huang, Y., Li, H. & Xiao, Y. Using 3dRPC for RNA–protein complex structure prediction. Biophys. Rep. 2, 95–99 (2016). 

45. Wild, K., Sinning, I. & Cusack, S. Crystal structure of an early proteinRNA assembly complex of the signal recognition particle. Science 294, 598–601 (2001). 

46. Palioura, S., Sherrer, R. L., Steitz, T. A., Söll, D. & Simonović, M. The human SepSecS-tRNASec complex reveals the mechanism of selenocysteine formation. Science 325, 321–325 (2009). 

47. Teplova,M. et al. Protein-RNA and protein-protein recognitionby dual KH1/2 domains of the neuronal splicing factor Nova-1. Structure 19, 930–944 (2011). 

48. Sekine, S. I. et al. ATP binding by glutamyl‐tRNA synthetase is switched to the productive mode by tRNA binding. EMBO J. 22, 676–688 (2003). 

49. Zhang, C., Shine, M., Pyle, A. M. & Zhang, Y. US-align: universal structure alignments of proteins, nucleic acids, and macromolecular complexes. Nat. Methods 19, 1109–1115 (2022). 

50. Liu, H., Jian, Y., Zeng, C. & Zhao, Y. RNA-protein interaction prediction using network-guided deep learning. https://doi.org/10. 5281/zenodo.14747845 (2025). 

## Acknowledgements 

This work was supported by the National Natural Science Foundation of China (grant no. 12175081); Hubei Science Fund for Distinguished Young Scholars (grant no. 2024AFA077); Fundamental Research Funds for the Central Universities (grant no. CCNU22QN004 and CCNU24JCPT011). The Central China Normal University’s excellent postgraduate education innovation funding project (grant no. 2024CXZZ146). 

## Author contributions 

H.L. built the RPI prediction algorithms, performed the analysis, and wrote the paper; Y.J. helped with the data collection; C.Z. helped with the network model discussion; Y.Z. supervised the overall study and edited the manuscript. 

## Competing interests 

The authors declare no competing interests. 

## Additional information 

Supplementary information The online version contains supplementary material available at https://doi.org/10.1038/s42003-025-07694-9. 

Correspondence and requests for materials should be addressed to Yunjie Zhao. 

Peerreviewinformation CommunicationsBiology thanksChandranNithin, Pranam Chatterjee and the other, anonymous, reviewer(s) for their contribution to the peer review of this work. Primary Handling Editors: Yuedong Yang and Laura Rodríguez Perez. 

Reprints and permissions information is available at http://www.nature.com/reprints 

Publisher’s note Springer Nature remains neutral with regard to jurisdictional claims in published maps and institutional affiliations. 

Open Access This article is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License, which permits any non-commercial use, sharing, distribution and reproduction in any medium or format, as long as you give appropriate credit to the original author(s) and the source, provide a link to the Creative Commons licence, and indicate if you modified the licensed material. You do not have permission under this licence to share adapted material derived from this article or parts of it. The images or other third party material in this article are included in the article’s Creative Commons licence, unless indicated otherwise in a credit line to the material. If material isnot includedin thearticle’s CreativeCommons licenceandyour intended use is not permitted by statutory regulation or exceeds the permitted use, you will need to obtain permission directly from the copyright holder. To view a copy of this licence, visit http://creativecommons.org/licenses/bync-nd/4.0/. 

© The Author(s) 2025 

Communications Biology | (2025) 8:247 

14 

