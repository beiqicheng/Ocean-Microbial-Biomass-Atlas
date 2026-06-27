# OMBA — Ocean Microbial Biomass Atlas

A modular Streamlit platform for microbial biomass estimation and validation.

## Structure

- Home / Project Manager
- SSU Explorer
- LSU Explorer
- RNA Yield
- Biomass
- Microscopy
- Quantitative Metagenome
- Environment
- Global Atlas

## File formats

`read_counts.txt`
- Sample
- SSU_hits
- LSU_hits

`SSU_matrix.txt` and `LSU_matrix.txt`
- Taxonomy
- Length
- Sample columns...

## Run locally
```bash
pip install -r requirements.txt
streamlit run app.py
```
