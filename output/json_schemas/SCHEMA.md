Schema of the Response of Json Response from SierraPy
=====================================================

This folder contains the json response from the sierrapy. This will help to parse the data well. Each json data that contains child structures is subdivided into further json schema so that we can get the outputs in a better way. the json schemas will be put in this file for better understanding

Overall JSON schemas
--------------------

```json
{
  "allGenes": [],
  "currentVersion": {},
  "currentProgramVersion": {},
  "report": {
    "inputSequence": {},
    "bestMatchingSubtype": {},
    "availableGenes": [],
    "mixtureRate": 0,
    "mutations": [],
    "unusualMutations": [],
    "frameShifts": [],
    "insertions": [],
    "deletions": [],
    "stopCodons": [],
    "ambiguousMutations": [],
    "apobecMutations": [],
    "mutationPrevalences": [],
    "algorithmComparison": [],
    "drugResistance": [],
    "alignedGeneSequences": [],
    "__typename": "SequenceAnalysis"
  }
}
```

However, here is the schema for this particular overall data.

```txt
allGenes: list
    []: list (empty)
currentVersion: dict
currentProgramVersion: dict
report: dict
    inputSequence: dict
    bestMatchingSubtype: dict
    availableGenes: list
        []: list (empty)
    mixtureRate: str
        int
    mutations: list
        []: list (empty)
    unusualMutations: list
        []: list (empty)
    frameShifts: list
        []: list (empty)
    insertions: list
        []: list (empty)
    deletions: list
        []: list (empty)
    stopCodons: list
        []: list (empty)
    ambiguousMutations: list
        []: list (empty)
    apobecMutations: list
        []: list (empty)
    mutationPrevalences: list
        []: list (empty)
    algorithmComparison: list
        []: list (empty)
    drugResistance: list
        []: list (empty)
    alignedGeneSequences: list
        []: list (empty)
    __typename: str
        str
```

Now let us go into details for the other parts of the json data.

-	Input sequence is just information about the sequence that we put in for analysis
-	Best matching subtype is the best subtype that the data came out to be. The stanford HIV DB runs analysis on the data and gives the particular subtype of this with the probability of its assurance.
-	The available genes talks about the genes that the algorithm finds in the nucleotide sequence that was put into it. It has the following structure.

	```txt
	name: str
	str
	__typename: str
	str
	```

A function that I suppose will be useful to get the available genes in this inputSequence

```python
availableGenes = [
{
"name":"gag",
"__typename":"Gene"
},
{
"name":"CA",
"__typename":"Gene"
},
{
"name":"pol",
"__typename":"Gene"
},
{
"name":"PR",
"__typename":"Gene"
},
{
"name":"RT",
"__typename":"Gene"
}
]

genes = [gene.name for gene in availableGenes]

# Output should be ['gag', 'CA', 'pol', 'PR', 'RT']
```

-	The mixtureRate must be something I don't know about the genomes.
-	The mutation aspect contains information on the gene and the mutation that occurs in the gene. It has multiple copies of this kind of information.

	```json
	{
	"gene":{
	"name":"PR",
	"__typename":"Gene"
	},
	"text":"I13V",
	"__typename":"Mutation"
	},
	```

We can access the necessary information from it by using the following line of code

```python
geneMutation = {dict(mutation.gene.name, mutation.text) for mutation in mutations }
```

-	Then we have the following. They are not so important but they can be present.

	```json
	"unusualMutations":[],
	"frameShifts":[],
	"insertions":[],
	"deletions":[],
	"stopCodons":[],
	"ambiguousMutations":[],
	"apobecMutations":[]
	```
