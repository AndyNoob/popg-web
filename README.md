# PopG Web (modified by AndyNoob)

The main algorithm is identical to [basecon's version](https://github.com/basecon/popg-web). However, the rendering of plot is disabled intentionally and the website now automatically generate a set of data for population sizes of 50, 100, 250, 500, 750, 1000, 1500, 2000, 2500, and 3000. Each set of data contains the number of generations before each population's allele within the set is fixed or lost. The data is outputted to the console when the website is loaded.
> [!NOTE]
> Each column represents a populate size. 
```
128	84	705	1249	1853	6248	3748	4086	2904	39659
57	81	92	714	2291	735	4980	1862	9595	7525
106	88	1063	2415	1678	896	3705	9518	4921	2569
66	214	399	1386	728	863	2759	15671	13880	4257
162	328	498	560	1407	8203	1879	3351	3150	1797
211	230	477	854	1716	2280	2580	3562	2576	13529
124	269	618	545	3074	3671	1752	30231	13894	20049
184	136	308	751	1611	2306	1243	2481	13485	5138
37	337	779	988	1533	3154	5322	24266	4489	14535
48	347	1403	358	433	1519	3355	17931	4680	2768
```