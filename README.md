# Translator Answer Appraiser
The Answer Appraiser is a micro-service that attaches a Confidence, Clinical Information, and Novelty score to each result of a given TRAPI message.

## Scoring Metrics
- The Confidence Score is currently a product of ARA scores that gives a strong weight to any ARA that scores a result very high.

- The Clinical Information Score is returned for queries where the results are drug-chemical pairs.  That factor takes into account the weighted average odds ratio of an observed association between a drug (query result) and a disease outcome (query) calculated across 3 independent knowledge providers.

- The Novelty Score captures many dimensions and holds the strong potential to differentiate results for researchers, beyond the use of a single and more conventional summarization based on published literature. These dimensions are contributed by a microservice termed the Annotator and include:
  - Publication recency.  There is a tradeoff in weighting highly cited earlier publications, versus newer publications that represent the latest novel findings, even if they have not had a chance to be cited widely yet. The recency factor evaluates the amount and year publications into a combined sigmoid function.  A recency factor near 1 denotes a query result that is supported by a very recent publication.
  - FDA Approval.  Drugs that are in a late phase (IV) of a US clinical trial are down weighted while those in earlier phases (I, II, III) get a greater novelty score.  Predicted drugs that are not in an FDA approval status are not scored and the approval status is specific to the disease being queried.
  - MolDist (Molecular Distinctiveness).  MolDist is calculated using the Jaccard coefficient between compounds shown in results and more established compounds.  Structurally identical compounds are scored as 1. 
  - TDL (Gene Target Development Levels).  Over 80% of the genes named in the results have a TDL: TClin, TChem, TBio or TDark, in decreasing levels of study knowledge of that gene (and thereby increasing novelty).
  - Gene distinctiveness. This metric denotes the number of pathways that a gene belongs to. 

  The Novelty score calculation starts with recency.  Then if the result is a chemical, MolSim and FDA approval may up or down-regulate the score.  Conversely, if the result is a gene, then the TDL and gene distinctiveness may up or down-regulate the score.


## Development Setup
Answer Appraiser can be run locally through either python or Docker.


### python
_Note: redis is also required and is easily stood up in docker._
1. Create and activate a virtual python environment (`python3.12 -m venv <path_to_venv>`, `source <path_to_venv>`)
1. `pip install -r requirements-lock.txt`
1. `python run.py`

### docker
1. `docker compose up`
