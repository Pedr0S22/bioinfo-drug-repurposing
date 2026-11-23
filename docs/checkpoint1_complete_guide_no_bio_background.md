# Checkpoint 1 — Complete Beginner-Friendly Guide (No Bioinformatics Background Assumed)

> **Purpose of this document:** Explain **every biological, chemical, statistical, and coding concept** behind our Checkpoint 1 work from scratch—using plain language, analogies, and our exact real numbers—so **every member of the group** can understand the entire pipeline and confidently defend it to the professor.

---

## Table of Contents

1. [The Big Picture: What Are We Actually Trying to Do?](#1-the-big-picture-what-are-we-actually-trying-to-do)
2. [Our Biological Target: What Is Acetylcholinesterase (AChE)?](#2-our-biological-target-what-is-acetylcholinesterase-ache)
3. [Why We Did NOT Use CDK2 or CDK9 (and How to Prove It)](#3-why-we-did-not-use-cdk2-or-cdk9-and-how-to-prove-it)
4. [The Databases We Used (and Why They Are 100% Open-Access)](#4-the-databases-we-used-and-why-they-are-100-open-access)
5. [Chemistry & Data Curation Basics: Turning Messy Lab Tables Into Clean Data](#5-chemistry--data-curation-basics-turning-messy-lab-tables-into-clean-data)
6. [How Molecules Become Numbers for Machine Learning](#6-how-molecules-become-numbers-for-machine-learning)
7. [Splitting Data Honestly: Random Split vs. Bemis–Murcko Scaffold Split](#7-splitting-data-honestly-random-split-vs-bemismurcko-scaffold-split)
8. [The #1 Trap in This Project: Data Leakage & Our 3-Tier Ablation](#8-the-1-trap-in-this-project-data-leakage--our-3-tier-ablation)
9. [All Checkpoint 1 Models & Metrics Explained (With Exact Calculations)](#9-all-checkpoint-1-models--metrics-explained-with-exact-calculations)
10. [Understanding Every Result We Got (Including Our Killer Error Analysis)](#10-understanding-every-result-we-got-including-our-killer-error-analysis)
11. [How to Explain the Two Plots & Notebooks in Class](#11-how-to-explain-the-two-plots--notebooks-in-class)
12. [How Our Code Reuses Assignment 5 & Works With `uv`](#12-how-our-code-reuses-assignment-5--works-with-uv)
13. [Professor Q&A Cheat Sheet (12 Tough Questions & Model Answers)](#13-professor-qa-cheat-sheet-12-tough-questions--model-answers)
14. [Mini-Glossary of Every Bioinformatics Term Used](#14-mini-glossary-of-every-bioinformatics-term-used)

---

## 1. The Big Picture: What Are We Actually Trying to Do?

### Drugs, Proteins, and the "Lock-and-Key" Idea
- Inside your body, **proteins** (especially **enzymes** and **receptors**) are tiny molecular machines that perform chemical jobs.
- Many diseases happen because a specific protein is **too active** (or not active enough).
- A **small-molecule drug** (like aspirin, caffeine, or ibuprofen—also called a **ligand**) is a tiny organic chemical designed to physically fit into a pocket on a specific protein—called the **target**—like a **key fitting into a lock**.
- When the drug fits into the pocket and blocks the protein from working, it is called an **inhibitor** (or more broadly, a **binder**).

### What Is Drug Repurposing?
- Developing a brand-new drug from scratch takes **10–15 years** and **$1–2 billion**, and ~90% fail because they turn out to be toxic or unsafe in humans.
- **Drug repurposing** asks a clever shortcut question: *"Among the ~2,000+ drugs that are **already FDA/EMA approved** (meaning we already know they are safe in humans and we already know how to manufacture them), could some of them secretly fit into the pocket of a **different** disease target?"*
- Testing 2,000 drugs in a wet lab is expensive and slow. **Bioinformatics / computational drug discovery** lets us score and **rank all ~2,000 approved drugs on a computer** so biologists only have to test the top 10–50 candidates in the lab.

### What Is a "Ligand-Based" Model (QSAR)?
There are two main ways to predict whether a drug binds a protein:
1. **Structure-based** (e.g., *molecular docking*): uses the **3D shape of the protein pocket** (from X-ray crystallography in the PDB) and tries to physically fit each 3D drug into the pocket.
2. **Ligand-based** (**QSAR** = *Quantitative Structure–Activity Relationship*): ignores the 3D protein completely! Instead, it looks only at **thousands of molecules previously tested in labs against that target** (some worked = *actives*, some didn't = *inactives*) and trains a machine learning model to recognize what chemical shapes/features make a molecule bind.
   - **Core assumption (Similarity Principle):** *Chemically similar molecules tend to have similar biological activity.*

### Project 7's Exact Logic (The "Sanity Test")
Anyone can train a model and print out a list of "top 10 repurposed drugs"—but how do you know the ranking isn't garbage?
- **The Positive-Control Test:** Among the ~2,000 approved drugs we rank, a small known subset is **already approved to treat our target** (our **positive controls**).
- **The Rule:** We **hide** all approved drugs from our model during training, then ask the model to rank all ~2,000 approved drugs from best to worst. **If the drugs already known to bind our target don't float near the very top of the ranking, our model cannot be trusted to propose new candidates.**

### Our One-Sentence Research Question
> *For human acetylcholinesterase (**AChE**, `CHEMBL220`), which FDA/EMA-approved small-molecule drugs do ligand-based QSAR models rank as candidate binders—and does a leak-free ranking recover the approved AChE drugs near the top of the library?*

---

## 2. Our Biological Target: What Is Acetylcholinesterase (AChE)?

You don't need a medical degree to explain AChE—just remember **three simple facts**:

1. **What AChE does in a healthy body:**
   - Nerve cells talk to each other (and to muscles) by releasing a chemical messenger called **acetylcholine (ACh)**.
   - Once the nerve signal is delivered, **Acetylcholinesterase (AChE)** is the enzyme that acts like a "broom," rapidly **breaking down acetylcholine** to turn the signal off.
2. **Why blocking (inhibiting) AChE is useful medically:**
   - In **Alzheimer's disease**, patients lose acetylcholine-producing neurons, so their memory/cognitive signals become too weak. Blocking AChE slows down the breakdown of acetylcholine, leaving more messenger available and improving symptoms. Famous approved Alzheimer's AChE drugs: **Donepezil**, **Rivastigmine**, **Galantamine**, **Tacrine**, **Benzgalantamine**.
   - In **Myasthenia gravis** (a neuromuscular weakness disease) and anesthesia reversal, blocking AChE strengthens nerve-to-muscle signals. Famous approved drugs: **Pyridostigmine**, **Neostigmine**, **Physostigmine**, **Edrophonium**, **Ambenonium**.
3. **Why AChE is a dream target for a bioinformatics project:**
   - It has **massive, clean experimental data** (`10,079` `IC50` records in ChEMBL).
   - It has **21 FDA/EMA-approved drugs** (`18` unique curated small molecules) spanning multiple chemical families—giving us strong statistical power for our positive-control test.
   - It has **76 high-resolution 3D crystal structures** in the PDB (e.g., `4M0E`, `6NTO`, `5HF5`), so in Week 3 we can even inspect its binding pocket or run optional docking on our top repurposing hypotheses!

---

## 3. Why We Did NOT Use CDK2 or CDK9 (and How to Prove It)

Initially, our group considered **CDK2** (`CHEMBL301`) and **CDK9** (`CHEMBL3116`). Here is why we checked them systematically and switched to **AChE**—and why this is actually a **huge strength** to present at Checkpoint 1!

### What Are Clinical Phases (`max_phase`)?
Before a drug can be prescribed to patients, it goes through clinical trials:
- **Phase 1:** Small safety test in humans.
- **Phase 2:** Does it seem to work in patients with the disease?
- **Phase 3:** Large randomized trial proving it works better than standard care.
- **Phase 4 (`max_phase = 4` / `APPROVAL`):** **Approved by FDA/EMA** and on the market!

### What Happened When We Checked CDK2 and CDK9 Live?
- **CDK2** and **CDK9** are *cyclin-dependent kinases* (proteins that control cell division and gene transcription in cancer).
- **Approved CDK drugs do exist—but they target CDK4 and CDK6** (*palbociclib, ribociclib, abemaciclib, trilaciclib* for breast cancer), **NOT CDK2 or CDK9**.
- Every drug specifically targeting CDK2 or CDK9 (*dinaciclib, alvocidib/flavopiridol, seliciclib, milciclib, tagtociclib*) is still stuck in **Phase 1, Phase 2, or Phase 3** clinical trials (many failed due to toxicity).
- Therefore, both **ChEMBL 37** and **Open Targets Platform** return **ZERO (`0`) Phase 4 approved drugs** whose primary mechanism of action is CDK2 or CDK9!
- **Why that breaks Project 7:** Imagine building a Google search engine for "best pizza places in Coimbra" and trying to score how well it works—except the official list of known pizza places is **empty (`0`)**. You literally cannot run Project 7's mandatory positive-control test (*"does the ranking recover the drugs already approved against that target?"*)!

### How We Selected AChE Objectively (13-Target Live Screen)
Instead of guessing, we wrote [scripts/screen_targets.py](../scripts/screen_targets.py) to query **every single target from Slide 13** plus **CDK2 and CDK9** live across ChEMBL 37, Open Targets Platform, and RCSB PDB:

| Target | Source List | ChEMBL ID | ChEMBL IC50 | ChEMBL Ki | Approved Drugs (ChEMBL ∩/∪ Open Targets) | 3D PDB Structures | Verdict |
|---|---|---|---:|---:|---:|---:|---|
| **AChE** | Slide 13 | `CHEMBL220` | **10,079** | 1,970 | **21** (16 ChEMBL, 18 OT) | **76** | **#1 Winner:** huge `IC50` dataset + 21 approved drugs + 76 PDBs |
| SLC6A4 | Slide 13 | `CHEMBL228` | 5,130 | 4,678 | 28 | 31 | Strong **Plan B1** |
| DRD2 | Slide 13 | `CHEMBL217` | 2,480 | 14,842 | 63 | 12 | Mostly `Ki` (not `IC50`) & few PDBs (**Plan B3**) |
| ABL1 | Slide 13 | `CHEMBL1862` | 5,295 | 1,781 | 8 | 71 | Strong **Plan B2** (kinase) |
| BRAF | Slide 13 | `CHEMBL5145` | 11,017 | 47 | 6 | 127 | Only 6 approved drugs |
| PARP1 | Slide 13 | `CHEMBL3105` | 4,993 | 1,212 | 6 | 102 | Only 6 approved drugs |
| ADORA2A | Slide 13 | `CHEMBL251` | 3,050 | 9,055 | 5 | 190 | Only 5 approved drugs; mostly `Ki` |
| NTRK1 / 2 / 3 | Slide 13 | `CHEMBL2815`… | 10,647 / 924 / 556 | 934 / 957 / 704 | 5 / 4 / 4 | 56 / 6 / 4 | Few approved drugs |
| DHFR | Slide 13 | `CHEMBL202` | 1,590 | 512 | 3 | 88 | Only 3 approved drugs (noisy metrics) |
| **CDK2** | Initial | `CHEMBL301` | 3,692 | 1,210 | **0** | 489 | **Ruled out:** 0 approved positive controls |
| **CDK9** | Initial | `CHEMBL3116` | 2,351 | 426 | **0** | 25 | **Ruled out:** 0 approved positive controls |

*(Bonus catch in the group's early draft `notebooks/codigo.ipynb`: 4 of its 6 hardcoded "CDK family" ChEMBL IDs were accidentally wrong targets—`CHEMBL302` was Rat Adenosine A2a, `CHEMBL359` was the bacterium *Bacillus subtilis*, `CHEMBL303` was Rat NMDA receptor, and `CHEMBL2972` was Human PFKP! Catching that early saved the whole project.)*

### How to Show the Professor Proof in 30 Seconds
- Open `https://platform.opentargets.org/target/ENSG00000123374` (CDK2) and `https://platform.opentargets.org/target/ENSG00000136807` (CDK9) → click **Drugs and Clinical Candidates** → zero Phase 4 / Approved entries.
- Open `https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL301/` (CDK2) and `CHEMBL3116` (CDK9) → **Drug Mechanisms** → zero Phase 4 entries.

---

## 4. The Databases We Used (and Why They Are 100% Open-Access)

Project rule: **"Use open-access data only."**

| Database | What It Is (Plain English) | What We Pull From It | Why Open-Access Compliant |
|---|---|---|---|
| **ChEMBL (Release 37)** | Europe's massive public database of lab experiments measuring how strongly chemicals bind proteins (run by EMBL-EBI). | 1. All human AChE (`CHEMBL220`) lab `IC50` measurements.<br>2. All approved small-molecule drugs (`max_phase = 4`).<br>3. Official approved mechanisms of action for AChE. | 100% open-access (CC BY-SA 3.0), free REST API, no login required. |
| **Open Targets Platform** | Public drug-target evidence platform combining ChEMBL, FDA, and clinical trials (explicitly cited on Slide 13!). | Cross-check of all approved (`APPROVAL` / `PHASE_4`) drugs targeting AChE (`ENSG00000087085`) so we don't miss any positive control. | 100% open-access GraphQL API, no login required. |
| **UniProt** | Global encyclopedia of protein sequences and biological functions. | Human AChE entry (`P22303`) for protein context and binding-site residues. | 100% open-access (CC BY 4.0). |
| **RCSB PDB** | Public archive of 3D atomic structures of proteins (often photographed with a drug sitting inside the pocket). | Count & top high-resolution crystal structures (`4M0E`, `6NTO`, `5HF5`) for optional Week 3 structural interpretation / docking. | 100% open-access (CC0). |

**Why did we use ChEMBL `max_phase = 4` instead of DrugBank?**  
Slide 13 mentions DrugBank, *however* DrugBank's full downloadable dataset now requires an academic registration and license agreement that restricts redistribution. Using **ChEMBL `max_phase = 4` small molecules + Open Targets** (explicitly recommended as Plan B in the group guide) guarantees **100% open-access compliance**, zero registration friction, and seamless matching with ChEMBL bioactivities.

---

## 5. Chemistry & Data Curation Basics: Turning Messy Lab Tables Into Clean Data

Raw biological databases are notoriously messy—different labs use different units, test the same drug multiple times, attach different salts, or make transcription errors. In Assignment 5, the professor called this: **"Curation is the modelling."**

### Concept A: What Is `IC50` and `pIC50`?
- **`IC50`** (*Half-Maximal Inhibitory Concentration*) = how much drug concentration you need to shut down **50%** of the enzyme's activity.
  - **Smaller `IC50` = MORE potent** (you need less drug to block the enzyme!).
  - Example: Drug A needs `10 nM` (nanomolar), Drug B needs `10,000 nM` (`10 µM`). Drug A is **1,000× stronger**!
- **`pIC50`** converts `IC50` onto a clean logarithmic scale where **BIGGER = MORE POTENT**:
  $$\text{pIC}_{50} = -\log_{10}(\text{IC}_{50}\text{ in mol/L}) = 9 - \log_{10}(\text{IC}_{50}\text{ in nM})$$

| `IC50` in Lab Units | `IC50` in `nM` | `pIC50` | Intuitive Meaning | Our Classification Label |
|---|---:|---:|---|---|
| `1 nM` | 1 | **9.0** | Ultra-potent (like a prescription drug) | **Active (`1`)** |
| `100 nM` | 100 | **7.0** | Very potent | **Active (`1`)** |
| `1 µM` | 1,000 | **6.0** | Solid hit / active cutoff | **Active (`1`)** (`pIC50 >= 6.0`) |
| `3 µM` | 3,000 | **5.5** | Borderline / noisy middle zone | **Dropped (grey zone)** (`5.0 < pIC50 < 6.0`) |
| `10 µM` | 10,000 | **5.0** | Weak / inactive cutoff | **Inactive (`0`)** (`pIC50 <= 5.0`) |
| `> 10 µM` (`>` censored) | > 10,000 | **< 5.0** | Lab reached max dose and STILL couldn't hit 50% inhibition! | **Rescued as Inactive (`0`)!** |

### Concept B: Why We Did NOT Mix `IC50` and `Ki`
- `IC50` depends slightly on experimental conditions (substrate concentration), whereas `Ki` is an intrinsic binding constant. Pooling them blindly into one `pIC50` column introduces systematic assay noise (Pitfall #5 in the group guide). Since AChE already has **9,882 binding `IC50` records**, we kept **binding assays (`assay_type == 'B'`) and `IC50` only**.

### Concept C: Why Censored Records (`> 10,000 nM`) Are Gold for Classification (and Dangerous for Regression)
- Sometimes a lab tests a molecule up to `10 µM` (`10,000 nM`) and it *never* blocks 50% of the enzyme, so they record `standard_relation = ">"` and `value = 10000 nM`.
- **In Regression** (predicting an exact continuous `pIC50` number): you **must drop** `>` records, because you don't know whether the true `pIC50` is `4.9`, `3.0`, or `0` (as you proved in Assignment 5 Q1c!).
- **In Binary Classification** (`Active >= 6.0` vs. `Inactive <= 5.0`): `IC50 > 10,000 nM` **mathematically guarantees** `pIC50 < 5.0`—so it is a **100% genuine experimental inactive!** Rescuing these gave us **`+585` extra real inactives** (fixing ChEMBL's bias toward publishing mostly active compounds!).

### Concept D: SMILES, Salts, Charge Neutralization, and `InChIKey14`
- **SMILES** (*Simplified Molecular-Input Line-Entry System*): a compact text string describing a molecule's atoms and bonds (e.g., `CCO` = ethanol, `c1ccccc1` = benzene ring).
- **Why strip salts (`largest_fragment_smiles`)?** Pills are often packaged as salts (e.g., *Donepezil hydrochloride* = `Donepezil.Cl`) so they dissolve better in your stomach. Once dissolved, the `.Cl` floats away—the **largest fragment** is the actual drug! If we don't strip salts, `Donepezil` and `Donepezil.Cl` look like two *different* strings, causing duplicate rows and **silent train/test data leakage**!
- **Why neutralize charges (`Uncharger`)?** The same molecule can be drawn protonated (`[NH3+]`) in one paper and neutral (`N`) in another. Neutralizing standardizes both to one canonical form.
- **What is `InChIKey14`?** An `InChIKey` is a standardized 27-character hash of a molecule's chemical structure (like a SHA-256 hash for chemistry). Its **first 14 characters** represent the exact **atom connectivity** (ignoring minor stereochemistry or protonation drawing quirks). Matching ChEMBL ↔ Approved Drugs on `InChIKey14` guarantees zero hidden duplicates!

### Concept E: Our Exact Dataset Check Funnel (What Went In, What Came Out)
1. **ChEMBL Bioactivities (`AChE`, `CHEMBL220`):**
   - Started with **`9,882`** raw binding `IC50` rows.
   - Removed **`629`** rows flagged as dubious by ChEMBL (`data_validity_comment`).
   - Removed **`1,993`** non-exact (`>`, `<`) rows, **`67`** unusable non-molar units (`%`, `ug/mL`), and **`1`** broken SMILES → **`7,192`** clean exact measurements.
   - Aggregated replicate measurements per molecule using the **median** → **`5,664`** unique exact-`pIC50` compounds (**Regression dataset**).
     - *Diagnostic (Assignment 5 Q1b):* `570` compounds were measured more than once, and **`47`** had replicates disagreeing by `> 2.0` log units (classic `µM` vs `nM` lab transcription errors—taking the **median** instead of mean/first row protects our labels!).
   - For **Classification**: `2,799` exact actives (`>= 6.0`) + `1,114` exact inactives (`<= 5.0`) + **`585`** rescued censored inactives (`> 10 µM`), dropping `1,751` grey-zone (`5–6`) compounds → after `InChIKey14` dedup, **`4,366`** unique compounds (`2,692` active = `61.7%`, `1,674` inactive = `38.3%`).
2. **Approved Drugs Library (to rank):**
   - Started with **`3,475`** raw ChEMBL `max_phase = 4` small-molecule records.
   - Removed invalid SMILES, inorganic salts/metals (e.g., lithium, cisplatin, contrast agents), tiny solvents (`MW < 100` or `< 6` heavy atoms), huge peptides/macrocycles (`MW > 900`), and salt duplicates (via `InChIKey14`).
   - Final curated library: **`2,027`** unique, drug-like organic approved small molecules, containing **`18`** matched positive controls (`3` of the `21` raw entries were salt/connectivity duplicates or outside the organic small-molecule filter).

---

## 6. How Molecules Become Numbers for Machine Learning

Machine learning models (like Logistic Regression or Random Forest) need **fixed-length numerical vectors**, not SMILES strings. We used the two standard cheminformatics representations from Assignment 5:

1. **ECFP4 Fingerprint (`2,048` binary bits — `0` or `1`):**
   - Think of it like **Bag-of-N-Grams / hashing trick for graphs**: around every atom in the molecule, RDKit looks at all neighbors up to **radius = 2 bonds** (diameter = 4 bonds → "ECFP**4**"), hashes each local chemical neighborhood into a slot from `0` to `2,047`, and sets that bit to `1`.
   - **Strength:** Captures exact local chemical motifs (rings, functional groups) that fit into a protein pocket.
   - **Weakness:** If a test molecule uses a completely new chemical skeleton, many of its bits were never seen in training!
2. **RDKit Physicochemical Descriptors (`12` continuous numbers):**
   - Global physical properties of the molecule:
     - `MolWt` (molecular weight / size)
     - `MolLogP` (lipophilicity = how greasy vs. water-loving it is)
     - `TPSA` (topological polar surface area = how much of the surface can form polar/H-bond interactions)
     - `NumHDonors`, `NumHAcceptors` (hydrogen-bond donors/acceptors—key "Velcro" contacts in pockets!)
     - `NumRotatableBonds` (flexibility), `HeavyAtoms`, `AromaticRings`, `FractionCSP3` (3D-ness vs flatness), `RingCount`, `QED` (overall drug-likeness 0–1), `LabuteASA` (surface area).
   - **Strength:** Smooth, transferable across different chemical families!
3. **ECFP4 + RDKit Descriptors (`2,060` features combined):**
   - Gives the model both fine-grained local substructures AND global physicochemical trends—and indeed beat either representation alone on our internal AChE validation!

### What Is Tanimoto Similarity?
When comparing two binary fingerprints $A$ and $B$ (each 2,048 bits), **Tanimoto similarity** is literally just the **Jaccard index**:
$$\text{Tanimoto}(A, B) = \frac{|A \cap B|}{|A \cup B|} = \frac{\text{bits that are 1 in BOTH molecules}}{\text{bits that are 1 in EITHER molecule}}$$
- `0.0` = share zero chemical substructures.
- `1.0` = identical ECFP4 fingerprint.
- `>= 0.85` = very close structural cousin ("analogue").

---

## 7. Splitting Data Honestly: Random Split vs. Bemis–Murcko Scaffold Split

Why can't you just use `train_test_split(..., shuffle=True)` in drug discovery?

### The Car Chassis Analogy
- Chemists don't invent 5,000 unrelated molecules—they pick one promising **core ring skeleton** (called the **Bemis–Murcko scaffold**, like a car chassis) and make 50–100 tiny variations by swapping side chains (changing the paint color or mirrors).
- **Random split** puts 35 variations of the *same* chassis in train and 15 in test. The model memorizes the chassis and looks like a genius (`ROC-AUC = 0.96`, `R² = 0.71`), even if it learned almost nothing that generalizes to new chemistry!
- **Bemis–Murcko Scaffold Split** (reused directly from Assignment 5 Task B) groups all molecules sharing the same core skeleton together and puts **whole scaffolds** either completely in train or completely in test—**0 scaffolds shared between train and test!**
- In our leak-free AChE dataset (`4,338` compounds), there are **`2,095` distinct Bemis–Murcko scaffolds** (largest series = `88` compounds, `1,492` singletons).
- Look how honest the scaffold split is: `RF Regressor (ECFP4+Desc)` drops from **`R² = 0.713` (random)** down to **`R² = 0.488` (scaffold)**, and `RF Classifier (ECFP4+Desc)` drops from **`ROC-AUC = 0.960` (random)** to **`0.900` (scaffold)**—still strong, but honest!

---

## 8. The #1 Trap in This Project: Data Leakage & Our 3-Tier Ablation

### What Is Data Leakage Here?
ChEMBL collects published lab measurements—so **of course** labs have measured the `IC50` of famous approved AChE drugs like Donepezil, Tacrine, and Galantamine!
- In our Dataset Check, we proved that:
  - **`14 out of 18` approved AChE drugs (78%!)** were sitting directly inside the raw ChEMBL regression set!
  - **`7 out of 18` approved AChE drugs (39%!)** were sitting directly inside the raw ChEMBL classification set!
- If you leave them in the training set, the model literally **memorizes the test answers** and trivially puts Donepezil and Tacrine at the top of the list—proving nothing!

### Our 3-Tier Leakage Ablation (A Killer Result for the Report)
Instead of only running one version, we ran **all three tiers** side by side to *measure* exactly how much leakage inflates performance:
1. **Tier 1 — Naive (Leaky):** Approved drugs left inside ChEMBL training (`5,664` reg / `4,366` clf).
2. **Tier 2 — Leak-Free (Exact) *(Our primary honest benchmark)*:** Removed every training compound matching **any** of the 2,027 approved drugs by standardized `InChIKey14` or canonical SMILES (`51` removed from reg → `5,613`; `28` removed from clf → `4,338`).
3. **Tier 3 — Strict Analogue-Free (`Tanimoto < 0.85`):** Also removed any training molecule that is a near-identical chemical cousin (`ECFP4 Tanimoto >= 0.85`) of any approved drug (`+3` removed from reg → `5,610`; `+1` from clf → `4,337`).

**Result:** Removing exact leakage dropped `RF Classifier` `EF@5%` from **`11.48 → 9.27`** (`BEDROC20: 0.546 → 0.426`) and `RF Regressor` `EF@5%` from **`8.61 → 7.07 → 6.18`** (`BEDROC20: 0.392 → 0.322 → 0.291`). Even after removing leakage, our models still dramatically beat random chance (`EF = 1.0`)!

---

## 9. All Checkpoint 1 Models & Metrics Explained (With Exact Calculations)

### The 7 Methods We Compared
We ran **both** problem formulations (classification and regression) plus a trivial sanity-check control, each averaged across **5 random seeds** (`0, 1, 2, 3, 4`) for honest error bars (`mean ± std`):

1. **Control: `MolWt alone` (Size-bias sanity check, Assignment 5 Q3b):** Ranks approved drugs solely by molecular weight (heaviest first). If this gets good enrichment, our problem is trivially biased by molecule size.
2. **Classification Baseline 1 — `Max Tanimoto (ECFP4)`:** Zero ML training! For each approved drug, compute its Tanimoto similarity to every *active* training compound and take the **maximum**. (Closest-cousin lookup—the classic cheminformatics baseline.)
3. **Classification Baseline 2 — `Logistic Regression (ECFP4)`:** Linear classifier (`class_weight="balanced"`) on the 2,048 binary ECFP4 bits; ranks approved drugs by predicted probability $P(\text{active})$.
4. **Classification Baseline 3 — `RF Classifier (ECFP4+Desc)`:** 150-tree Random Forest classifier (`max_features="sqrt"`, `class_weight="balanced"`) on ECFP4 + 12 RDKit descriptors.
5. **Regression Baseline 1 — `RF Regressor (ECFP4)`:** 150-tree Random Forest predicting continuous `pIC50` from ECFP4; ranks approved drugs by highest predicted `pIC50`.
6. **Regression Baseline 2 — `RF Regressor (RDKit Desc)`:** Same, using only the 12 physicochemical descriptors.
7. **Regression Baseline 3 — `RF Regressor (ECFP4+Desc)`:** Same, combining both representations (exact Assignment 5 Part 3 setup!).

### Every Ranking Metric Explained Simply (Using Our Real Numbers!)
We are ranking **`N = 2,027` approved drugs**, of which **`18` are positive controls** (`18 / 2,027 = 0.89%` base rate—very imbalanced!).

1. **ROC-AUC** *(0.50 = random, 1.00 = perfect)*:
   - Probability that a randomly chosen approved AChE drug gets a higher score than a randomly chosen non-AChE drug.
   - *Watch out:* Because 99.1% of drugs are negatives, ROC-AUC gives equal weight to moving a drug from rank #1,500 to #800 as moving it from #50 to #5—even though biologists only ever test the **very top** of the list!
2. **PR-AUC** *(Precision-Recall AUC; random = base rate = **`0.0089`**, 1.00 = perfect)*:
   - Specifically measures whether top-scoring predictions are genuinely positive under severe class imbalance.
   - Our leak-free `Max Tanimoto` gets **`PR-AUC = 0.467`** (**52× higher than random `0.0089`**!), and `RF Classifier` gets **`0.247`** (**28× random**).
3. **Enrichment Factor (`EF@1%` and `EF@5%`)** *(reused from Assignment 5 Task C; `1.0` = random)*:
   - **"How many times more real AChE drugs are in the top X% of my list than random picking would put there?"**
   - Formula: $\text{EF@}f = \frac{\text{hits in top } k / k}{\text{total positives} / N}$, where $k = \lceil f \times N \rceil$.
   - **Concrete example using our real `Max Tanimoto` leak-free numbers:**
     - Top `1%` of `2,027` drugs = $k = \lceil 0.01 \times 2027 \rceil = \mathbf{21\text{ drugs}}$.
     - `Max Tanimoto` put **`9` of the `18` real AChE drugs** inside the top 21!!
     - Hit rate in top 1% = $9 / 21 = 42.9\%$. Random hit rate = $18 / 2027 = 0.888\%$.
     - $\text{EF@1\%} = \frac{9 / 21}{18 / 2027} = \mathbf{48.26\times}$ better than random!
     - Top `5%` = $k = \lceil 0.05 \times 2027 \rceil = \mathbf{102\text{ drugs}}$, containing **`12 / 18` (`66.7%`)** real AChE drugs → $\text{EF@5\%} = \mathbf{13.25\times}$ better than random!
4. **BEDROC ($\alpha = 20$)** *(0.0 = worst, ~0.05 = random, 1.0 = perfect)*:
   - *Boltzmann-Enhanced Discrimination of ROC*. Like ROC-AUC, bounded cleanly between `0` and `1`, **except** it uses an exponential decay weight ($\alpha = 20$) so **80% of the score comes from the top 8% of the ranked list**! This is the gold-standard "early recognition" metric in virtual screening.
5. **Median Rank & % in Top 5%**:
   - The easiest number to explain out loud: *"Out of 2,027 approved drugs, half of all approved AChE drugs appeared in the top **21** (`Max Tanimoto`) or top **155** (`RF Classifier`), compared to rank **1,014** by random chance!"*
6. **Applicability Domain (`AD_max_tanimoto`)**:
   - For every approved drug, we compute its maximum Tanimoto similarity to our leak-free training actives.
   - High (`> 0.5–0.8`) = **in-domain** (chemically related to known AChE binders; predictions are trustworthy).
   - Low (`< 0.3`) = **out-of-domain** (chemically alien to training data; ML predictions may be risky extrapolations—OR exciting *scaffold hops* if supported by independent evidence!).

---

## 10. Understanding Every Result We Got (Including Our Killer Error Analysis)

### Summary Table — Leak-Free (Exact) Approved-Drug Ranking (`N = 2,027`, `18` Positive Controls)

| Formulation | Method | ROC-AUC | PR-AUC | EF@1% | EF@5% | BEDROC20 | Median Rank | % Top 5% |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Control | **MolWt alone** | 0.378 | 0.007 | **0.00** | **0.00** | 0.018 | 1,440 | 0.0% |
| Classification | **Max Tanimoto (ECFP4)** | **0.847** | **0.467** | **48.26** | **13.25** | **0.654** | **21.0** | **66.7%** |
| Classification | Logistic Regression (ECFP4) | 0.654 | 0.021 | 5.36 | 3.31 | 0.149 | 635.0 | 16.7% |
| Classification | **RF Classifier (ECFP4+Desc)** | 0.631 | 0.247 | 27.88 | **9.27** | 0.426 | 154.5 | **46.7%** |
| Regression | RF Regressor (ECFP4) | 0.540 | 0.133 | 16.09 | 3.53 | 0.198 | 812.8 | 17.8% |
| Regression | RF Regressor (RDKit Desc) | 0.589 | 0.032 | 5.36 | 4.42 | 0.178 | 700.2 | 22.2% |
| Regression | **RF Regressor (ECFP4+Desc)** | 0.613 | 0.200 | 19.30 | **7.07** | 0.322 | 522.9 | **35.6%** |

### 5 Big Scientific Insights to Highlight in Your Presentation

1. **Size Bias Is Completely Ruled Out:**
   - `MolWt alone` gets `EF@1% = 0.00`, `EF@5% = 0.00`, `ROC-AUC = 0.378` (worse than random—AChE drugs tend to be slightly smaller than average approved drugs!). So our models are genuinely learning chemical binding motifs, not molecule size.
2. **Why Did Simple `Max Tanimoto` Beat ML Models in Checkpoint 1?**
   - Even after removing exact approved drugs, ChEMBL contains rich medicinal-chemistry series derived from classic AChE drug scaffolds (tacrine, donepezil, galantamine, carbamates). So `Max Tanimoto` directly recovered **Tacrine (#1), Donepezil (#2), Physostigmine (#3), Galantamine (#5), Rivastigmine (#6), Ambenonium (#7), Benzgalantamine (#9), Edrophonium (#13), Minaprine (#20), Neostigmine (#22)**!
   - **Why is this a GREAT Checkpoint 1 result?** It sets an **honest, tough baseline** (`EF@5% = 13.25`, `BEDROC20 = 0.654`) that our Week 2 challenger models (tuned XGBoost, Tanimoto-kernel SVM, ChemBERTa embeddings, calibrated ensembles) must try to beat—especially on **novel scaffolds**!
3. **Where ML Already Beats `Max Tanimoto` (The Scaffold-Hop Case!):**
   - Look at **Hexafluorenium** (an approved AChE inhibitor with low structural similarity to ChEMBL training actives: `AD_max_tanimoto = 0.27`). `Max Tanimoto` completely missed it (**rank `#802`**), whereas `RF Regressor (ECFP4+Desc)` rescued it to **rank `#32`** (top 1.6%!) by using transferable physicochemical descriptors + subtree features!
4. **Why Classification Beat Regression on Approved-Drug Enrichment (`EF@5%: 9.27 vs. 7.07`):**
   - Rescuing **`+585` censored experimental inactives (`IC50 > 10 µM`)** gave the classifier real negative boundary examples that exact-value regression had to discard.
   - Meanwhile, uncalibrated **Logistic Regression** struggled in the top 1% (`EF@1% = 5.36`) because linear weights on 2,048 sparse bits can assign extreme probabilities to unrelated out-of-domain drugs (its top 10 false positives all had `AD_max_tanimoto ≈ 0.23–0.31`).
5. **The "Showstopper" Biological Error Analysis (Why 3 Approved AChE Drugs Ranked Near the Bottom in EVERY Model):**
   - **Isoflurophate** (rank ~1,996, `AD = 0.14`) and **Echothiophate** (rank ~1,297, `AD = 0.23`): These are **irreversible covalent organophosphate inhibitors** (tiny phosphorus-containing eye drops related to nerve-agent chemistry that chemically react and form a permanent bond with AChE). Almost all ChEMBL AChE medicinal chemistry consists of *reversible non-covalent* inhibitors—so no ligand-based model trained on reversible drugs could ever recognize covalent organophosphates!
   - **Pralidoxime** (rank ~1,016–2,026, `AD = 0.25`): **Pralidoxime is NOT an AChE inhibitor at all!** It is an **AChE reactivator**—an emergency antidote that *repairs* AChE after organophosphate poisoning! ChEMBL/Open Targets automatically tag it as "acting on AChE," so it got swept into our positive controls. Discovering this through error analysis—and proposing to filter positive controls by `action_type == 'INHIBITOR'` in Week 2—shows true scientific rigor!

---

## 11. How to Explain the Two Plots & Notebooks in Class

### Plot 1: [ache_dataset_check.png](../results/figures/ache_dataset_check.png) *(from Notebook 01)*
- **Left panel (Histogram of curated exact `pIC50`):** Shows the distribution of all `5,664` unique compounds. The orange dashed line (`5.0`) and green dashed line (`6.0`) show our inactive/active cutoffs, with the `1,751` ambiguous grey-zone compounds sitting between them.
- **Middle panel (Stacked bar chart of Leak-Free Classification Balance):** Shows `2,692` actives vs. `1,674` inactives—highlighting in orange the **`+585` rescued censored inactives (`> 10 µM`)** that meaningfully improved class balance.
- **Right panel (Top 20 Bemis–Murcko Scaffold Series):** Shows the classic "long-tail" of medicinal chemistry: a few heavily explored scaffolds (largest = 88 analogues) followed by `1,492` singletons—visually proving why a **scaffold split** is mandatory!

### Plot 2: [ache_checkpoint1_baselines.png](../results/figures/ache_checkpoint1_baselines.png) *(from Notebook 02)*
- **Left panel (Grouped bar chart — 3-Tier Leakage Ablation on `EF@5%`):** Shows orange (*Naive/Leaky*) > blue (*Leak-Free Exact*) ≥ green (*Strict `Tanimoto < 0.85`*), visually quantifying how much leakage inflates performance, while `MolWt` sits flat at `0.0` (below the dashed random line at `1.0`).
- **Right panel (Leak-Free Cumulative Enrichment Curve, top 0–30% of library):** X-axis = % of the 2,027 approved drugs screened; Y-axis = % of the 18 real AChE drugs recovered. The steep jump of `Max Tanimoto` (orange), `RF Classifier` (green), and `RF Regressor` (purple) straight up the left edge shows dramatic early enrichment compared to the grey diagonal (*Random*) and `MolWt` (dark grey).

---

## 12. How Our Code Reuses Assignment 5 & Works With `uv`

We organized the repository cleanly so nothing is copy-pasted spaghetti:
- **`pyproject.toml` + `uv.lock`:** Run `uv sync` once and everyone gets the **exact same Python 3.12 environment** (`.venv`) with RDKit, scikit-learn, pandas, matplotlib, seaborn, and Jupyter.
- **`src/bioinfo_repurposing/`** directly wraps and extends the functions we passed in **Assignment 5 (CAP5)**:
  - `curation.py` ← Assignment 5 **Task A** (`largest_fragment_smiles`, `curate_activity`, replicate-disagreement check) + `Uncharger`, `InChIKey14`, censored inactive rescue, and `remove_leakage`.
  - `features.py` ← Assignment 5 **Task B** (`murcko_scaffold`, `scaffold_split`, `random_split`, `ecfp4`, `descriptors`) + `max_tanimoto_similarity`.
  - `metrics.py` ← Assignment 5 **Task C** (`enrichment_factor`) + `bedroc_score`, `ROC-AUC`, `PR-AUC`, and rank summaries.
  - `data_fetch.py` ← Clean, cached REST/GraphQL clients for ChEMBL 37, Open Targets v4, and RCSB PDB.
  - `baselines.py` ← Reproducible 5-seed classification + regression + 3-tier leakage pipeline.
- **`tests/test_core.py`:** Automated unit tests proving our `curate_activity`, `scaffold_split`, and `enrichment_factor` pass Assignment 5's exact assertions.

---

## 13. Professor Q&A Cheat Sheet (12 Tough Questions & Model Answers)

1. **"Wait—didn't you originally pick CDK2 or CDK9? Why did you change to AChE?"**  
   *Answer:* "We ran a live dataset check across ChEMBL 37 and Open Targets on CDK2, CDK9, and all Slide 13 targets. Both databases confirmed **0 FDA-approved Phase 4 drugs** with primary mechanism of action on CDK2 or CDK9 (approved CDK drugs target CDK4/6). Because Project 7's core validation requires recovering already-approved drugs for the target, we switched to **AChE (`CHEMBL220`)**, which has 10,079 single-endpoint `IC50` records, 21 approved positive controls, and 76 holo PDB structures."
2. **"How do you know you didn't leak approved drugs into your training set?"**  
   *Answer:* "We didn't match only by ChEMBL ID or raw SMILES—we stripped salts, neutralized charges, and matched on the 14-character **InChIKey connectivity layer**. That caught **14/18 positive controls (51 approved drugs total)** in the raw regression set and **7/18 (28 total)** in classification. We also ran a stricter `Tanimoto < 0.85` analogue-removal tier and reported all three tiers side by side."
3. **"Why did you choose `pIC50 >= 6.0` as active and `<= 5.0` as inactive?"**  
   *Answer:* "`pIC50 = 6.0` (`1 µM`) is the standard medicinal-chemistry hit potency cutoff, and `<= 5.0` (`>= 10 µM`) clearly separates weak/inactive compounds. Dropping the `5–6` grey zone prevents borderline assay noise (~0.3–0.5 log units) from flipping class labels—and we *also* ran continuous `pIC50` regression using all `5,664` exact compounds so we could compare both formulations directly."
4. **"Why did you keep `>` censored measurements in classification when Assignment 5 dropped them?"**  
   *Answer:* "Assignment 5 did *regression*, where encoding `> 10 µM` as `pIC50 = 5.0` biases the squared error loss. In *binary classification* with an inactive threshold of `<= 5.0`, `IC50 > 10,000 nM` is mathematically guaranteed to be `pIC50 < 5.0`. Rescuing those `+585` compounds gave us real experimental negatives instead of synthetic decoys."
5. **"Why did you drop `Ki` records instead of combining `IC50` and `Ki`?"**  
   *Answer:* "`IC50` is assay-condition dependent whereas `Ki` is an intrinsic dissociation constant, so mixing them without a Cheng–Prusoff correction adds noise. Because AChE already has nearly 10,000 binding `IC50` records, sticking strictly to `IC50` keeps the label physically consistent."
6. **"Why did you use ChEMBL `max_phase = 4` instead of DrugBank?"**  
   *Answer:* "DrugBank's full download requires registration and an academic license restricting redistribution, which could conflict with the open-access-only rule. ChEMBL `max_phase = 4` + Open Targets Platform are 100% open-access, cited on Slide 13, and reproducible via public APIs."
7. **"Why is your scaffold-split performance lower than your random-split performance?"**  
   *Answer:* "Random splits put close analogues from the same Bemis–Murcko scaffold on both sides of the split, allowing models to memorize chemical series. Our scaffold split holds out whole chemical families (0 shared scaffolds), giving an honest estimate of generalization to novel chemistry."
8. **"Why report `EF@1%`, `EF@5%`, `PR-AUC`, and `BEDROC20` instead of just `ROC-AUC` or accuracy?"**  
   *Answer:* "Only 18 out of 2,027 approved drugs (0.89%) are positive controls, so accuracy would be 99.1% for a useless all-zero model, and ROC-AUC weights the whole list equally. Early-enrichment metrics (`EF@1%/5%`, `BEDROC20`, `PR-AUC`) specifically measure whether true binders appear at the very top where a lab would actually test them."
9. **"Your simplest baseline (`Max Tanimoto`) beat Random Forest and Logistic Regression—why?"**  
   *Answer:* "Even after exact deduplication, ChEMBL contains rich reversible AChE series structurally related to drugs like donepezil, tacrine, and galantamine, so nearest-neighbor similarity excels in-domain (`EF@5% = 13.25`). However, `Max Tanimoto` fails on low-similarity scaffold hops—for example, it ranked **Hexafluorenium** `#802` (`AD = 0.27`), whereas `RF Regressor` ranked it `#32`. In Week 2 we will test whether tuned XGBoost, Tanimoto-kernel SVM, or ChemBERTa embeddings can match Tanimoto in-domain *and* beat it on scaffold hops."
10. **"Why did Isoflurophate, Echothiophate, and Pralidoxime rank near the bottom?"**  
    *Answer:* "Isoflurophate and Echothiophate are irreversible *covalent* organophosphates chemically unrelated to reversible ChEMBL inhibitors, and **Pralidoxime** is actually an AChE *reactivator* antidote (not an inhibitor!). Filtering positive controls to true reversible inhibitors is one of our immediate Week 2 refinements."
11. **"Could your high enrichment just be molecular-weight bias?"**  
    *Answer:* "No—we explicitly included a `MolWt alone` control (as in Assignment 5), which scored `EF@1% = 0.00`, `EF@5% = 0.00`, and `ROC-AUC = 0.378`."
12. **"What is your Plan B if another group already chose AChE?"**  
    *Answer:* "Because our pipeline is completely parameterized by target symbol, we can switch in minutes to **`SLC6A4`** (5,130 IC50, 28 approved drugs, 31 PDBs) or **`ABL1`** (5,295 IC50, 8 approved drugs, 71 PDBs), both already screened in Notebook 01."

---

## 14. Mini-Glossary of Every Bioinformatics Term Used

| Term | Plain-English Meaning |
|---|---|
| **Target** | The human protein (enzyme/receptor) a drug binds to produce a medical effect (ours: **AChE**). |
| **Ligand** | A small molecule (drug or experimental compound) that binds to a target protein. |
| **Inhibitor vs. Reactivator** | An inhibitor blocks an enzyme; a reactivator (like Pralidoxime) restores a chemically poisoned enzyme. |
| **Bioactivity / `IC50` / `pIC50`** | Measured binding/inhibition strength. `pIC50 = 9 - log10(IC50 in nM)`; higher = stronger (`>= 6` = active). |
| **Censored measurement (`>`)** | Lab hit its max test concentration without seeing 50% inhibition (e.g., `> 10,000 nM` = guaranteed weak/inactive). |
| **SMILES** | Text notation representing a molecule's 2D graph of atoms and bonds. |
| **Salt / Counter-ion** | Extra ion packaged with a drug pill (e.g., `.HCl`); stripped during curation so only the active parent drug remains. |
| **InChIKey (14-char)** | Standardized chemical hash whose first 14 characters encode exact atom connectivity (used to prevent leakage). |
| **ECFP4 (Morgan fingerprint)** | 2,048-bit binary vector summarizing all local circular atom neighborhoods up to 2 bonds away. |
| **RDKit Descriptors** | Continuous physicochemical numbers describing a whole molecule (weight, logP/greasiness, polarity, rings, QED). |
| **Tanimoto Similarity** | Jaccard overlap (`0` to `1`) between two molecular fingerprints. |
| **Bemis–Murcko Scaffold** | Core ring-and-linker skeleton of a molecule after stripping side chains (used for honest chemical splits). |
| **Applicability Domain (AD)** | Region of chemical space close enough to training data for model predictions to be trustworthy. |
| **Positive Control** | Drugs *already approved* for our target, used to test whether our ranking puts known truths near the top. |
| **Enrichment Factor (`EF@X%`)** | How many times more positive controls appear in the top `X%` of our ranking than random guessing (`1.0`). |
| **BEDROC ($\alpha = 20$)** | Early-recognition ranking score (`0` to `1`) emphasizing the top ~8% of the ranked library. |
| **Holo PDB Structure** | Experimental 3D structure of a protein with a ligand bound inside its pocket (useful for optional Week 3 docking). |
