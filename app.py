import streamlit as st

st.set_page_config(page_title="OMBA", layout="wide")
st.title("Ocean Microbial Biomass Atlas (OMBA)")
st.caption("Project Manager and entry point for OMBA.")

st.info("This scaffold is ready for the first pages: SSU Explorer, LSU Explorer, RNA Yield, Biomass, Microscopy, Quantitative Metagenome, and Environment.")

st.subheader("Project Manager")
st.text_input("Project name", value=st.session_state.get("project_name", "AMT30"), key="project_name")

st.write("Upload your project files once and reuse them across pages.")
st.write("- read_counts.txt")
st.write("- SSU_matrix.txt")
st.write("- LSU_matrix.txt")
st.write("- optional RNA, microscopy, and environment tables")

st.success("Next: add the taxonomy browser and save the parsed tables under the project name.")
