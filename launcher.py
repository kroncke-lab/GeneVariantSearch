import streamlit as st

st.set_page_config(
    page_title="PubMed Tools Suite",
    page_icon="🔬",
    layout="centered"
)

st.title("🔬 PubMed Genetic Variant Tools")
st.markdown("---")

col1, col2 = st.columns(2)

with col1:
    st.markdown("### 🔍 Variant Extractor")
    st.markdown("""
    **Search PubMed** for genetic variants and extract clinical data from:
    - Abstracts
    - Full-text articles
    - Supplemental files (PDF, Excel, Word)
    
    **Use for:** Discovery and data extraction
    """)
    
    if st.button("Launch Extractor →", use_container_width=True, type="primary"):
        st.info("Run this command in Shell:")
        st.code("streamlit run app.py --server.port 5000", language="bash")
        st.markdown("Then refresh this page to see the app")

with col2:
    st.markdown("### ✅ Variant Validator")
    st.markdown("""
    **Upload a CSV** with annotated variants and validate against source papers.
    
    Checks if your annotations match what's actually reported.
    
    **Use for:** Quality control and verification
    """)
    
    if st.button("Launch Validator →", use_container_width=True, type="primary"):
        st.info("Run this command in Shell:")
        st.code("streamlit run validator_app/app.py --server.port 5000", language="bash")
        st.markdown("Then refresh this page to see the app")

st.markdown("---")
st.markdown("""
### 📝 Quick Start

**For Variant Extraction:**
```bash
streamlit run app.py --server.port 5000
```

**For Variant Validation:**
```bash
streamlit run validator_app/app.py --server.port 5000
```

**Switching between apps:** Stop the current workflow and run the other command.
""")
