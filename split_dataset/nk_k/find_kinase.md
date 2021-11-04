1. pdbid ---> uniport id ---> EC number. 
    All pdbid in PDBbind v2019 general set is converted
    According to https://en.wikipedia.org/wiki/List_of_EC_numbers_(EC_2) 
    EC number 2.7.1 --- 2.7.4, 2.7.9---2.7.14 and 2.7.99 are considered as kinase

2. KLIFS
    Kinase download from KLIF database

3. Human_search
    If the protein name containing "kinase" in PDBbind dataset is considered as kinase

4. Human_check
    If the pdbid is shared among above three source, it is collected as kinase, other pdbid is checked carefully.

5. Results
|subset|kinase|nonkinase|
|---|---|---|
|Refine|491|4361|
|General|
