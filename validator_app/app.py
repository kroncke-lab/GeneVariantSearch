import streamlit as st
import pandas as pd
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pubmed_api import PubMedAPI
from validator import VariantValidator

st.set_page_config(
    page_title="PubMed Variant Validator",
    page_icon="✅",
    layout="wide"
)

st.title("✅ PubMed Variant Annotation Validator")
st.markdown("""
Upload a CSV with annotated genetic variants and PMIDs. This tool will:
1. Fetch each paper from PubMed
2. Extract what the paper actually reports about the variant
3. Compare against your annotations
4. Output: **Correct** ✓ / **Incorrect** ✗ / **Unsure** ? with explanations
""")

with st.sidebar:
    st.header("⚙️ Configuration")
    
    model_choice = st.selectbox(
        "Select AI Model",
        ["Google Gemini 2.5 Flash (Free)", "Anthropic Claude 3 Haiku", "OpenAI GPT-4o Mini"],
        help="Gemini offers free tier at https://aistudio.google.com/apikey"
    )
    
    api_key = None
    if "Gemini" in model_choice:
        api_key = os.getenv("GEMINI_API_KEY")
        if not api_key:
            api_key = st.text_input("Gemini API Key", type="password", help="Get free key at https://aistudio.google.com/apikey")
    elif "Claude" in model_choice:
        api_key = os.getenv("ANTHROPIC_API_KEY")
        if not api_key:
            api_key = st.text_input("Anthropic API Key", type="password")
    else:
        api_key = os.getenv("OPENAI_API_KEY")
        if not api_key:
            api_key = st.text_input("OpenAI API Key", type="password")
    
    st.divider()
    st.markdown("### 📋 CSV Format")
    st.markdown("""
    Your CSV should contain:
    - **PMID** column (required)
    - **Variant** column (HGVS or protein name)
    - **Phenotype/Diagnosis** columns
    
    The app will auto-detect column names.
    """)

uploaded_file = st.file_uploader("Upload CSV with variant annotations", type=['csv'])

if uploaded_file:
    try:
        df = pd.read_csv(uploaded_file)
        
        st.success(f"✓ Loaded {len(df)} rows with {len(df.columns)} columns")
        
        with st.expander("📊 Preview uploaded data"):
            st.dataframe(df.head(10))
        
        st.subheader("🔍 Column Mapping")
        st.markdown("Help us identify the key columns in your CSV:")
        
        col1, col2 = st.columns(2)
        
        with col1:
            pmid_col = st.selectbox(
                "PMID Column",
                options=df.columns.tolist(),
                index=next((i for i, col in enumerate(df.columns) if 'pmid' in col.lower()), 0)
            )
        
        with col2:
            gene_col = st.selectbox(
                "Gene Column (optional)",
                options=["None"] + df.columns.tolist(),
                index=next((i+1 for i, col in enumerate(df.columns) if 'gene' in col.lower()), 0)
            )
        
        col3, col4 = st.columns(2)
        
        with col3:
            variant_cols = st.multiselect(
                "Variant Column(s)",
                options=df.columns.tolist(),
                default=[col for col in df.columns if any(keyword in col.lower() for keyword in ['variant', 'hgvs', 'protein', 'mutation'])]
            )
        
        with col4:
            phenotype_cols = st.multiselect(
                "Phenotype/Diagnosis Column(s)",
                options=df.columns.tolist(),
                default=[col for col in df.columns if any(keyword in col.lower() for keyword in ['phenotype', 'cancer', 'diagnosis', 'diag', 'lqt', 'disease'])]
            )
        
        if not api_key:
            st.warning("⚠️ Please add your API key in the sidebar to continue")
        elif not pmid_col:
            st.error("❌ Please select a PMID column")
        elif not variant_cols:
            st.error("❌ Please select at least one variant column")
        else:
            max_rows = st.slider("Maximum rows to validate", min_value=1, max_value=min(100, len(df)), value=min(10, len(df)))
            
            if st.button("🚀 Start Validation", type="primary"):
                validator = VariantValidator(api_key, model_choice)
                
                results = []
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                df_subset = df.head(max_rows)
                
                for idx, row in df_subset.iterrows():
                    row_num = int(idx) + 1
                    progress = row_num / max_rows
                    progress_bar.progress(progress)
                    
                    pmid = str(row[pmid_col]).strip()
                    if not pmid or pmid.lower() == 'nan':
                        status_text.text(f"Row {row_num}/{max_rows}: Skipping (no PMID)")
                        results.append({
                            'row_index': int(idx),
                            'pmid': pmid,
                            'validation_status': 'SKIPPED',
                            'validation_reason': 'No PMID provided',
                            'paper_title': '',
                            'reported_data': ''
                        })
                        continue
                    
                    gene_info = str(row[gene_col]) if gene_col != "None" and gene_col in df.columns and pd.notna(row[gene_col]) else ""
                    variant_info = " | ".join([str(row[col]) for col in variant_cols if pd.notna(row[col])])
                    phenotype_info = " | ".join([str(row[col]) for col in phenotype_cols if pd.notna(row[col])])
                    
                    status_text.text(f"Row {row_num}/{max_rows}: Validating PMID {pmid}...")
                    
                    result = validator.validate_variant(
                        pmid=pmid,
                        gene=gene_info,
                        annotated_variant=variant_info,
                        annotated_phenotype=phenotype_info
                    )
                    
                    result['row_index'] = int(idx)
                    results.append(result)
                
                progress_bar.progress(1.0)
                status_text.text("✓ Validation complete!")
                
                results_df = pd.DataFrame(results)
                output_df = df_subset.copy()
                output_df['validation_status'] = results_df['validation_status']
                output_df['validation_reason'] = results_df['validation_reason']
                output_df['paper_title'] = results_df['paper_title']
                output_df['reported_data'] = results_df['reported_data']
                
                st.subheader("📊 Validation Results")
                
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    correct_count = (output_df['validation_status'] == 'CORRECT').sum()
                    st.metric("✅ Correct", correct_count)
                with col2:
                    incorrect_count = (output_df['validation_status'] == 'INCORRECT').sum()
                    st.metric("❌ Incorrect", incorrect_count)
                with col3:
                    unsure_count = (output_df['validation_status'] == 'UNSURE').sum()
                    st.metric("❓ Unsure", unsure_count)
                with col4:
                    skipped_count = (output_df['validation_status'] == 'SKIPPED').sum()
                    st.metric("⊘ Skipped", skipped_count)
                
                st.dataframe(
                    output_df,
                    column_config={
                        'validation_status': st.column_config.TextColumn(
                            'Status',
                            help='Validation result'
                        ),
                        'validation_reason': st.column_config.TextColumn(
                            'Reason',
                            help='Explanation for validation result'
                        )
                    }
                )
                
                csv_output = output_df.to_csv(index=False)
                st.download_button(
                    label="📥 Download Validated CSV",
                    data=csv_output,
                    file_name="validated_variants.csv",
                    mime="text/csv"
                )
                
    except Exception as e:
        st.error(f"Error loading CSV: {str(e)}")
        st.exception(e)
