import streamlit as st
import pandas as pd
from pubmed_api import PubMedAPI
from llm_analyzer import LLMAnalyzer
import os
from datetime import datetime
import sys

# Check if running through streamlit CLI
try:
    from streamlit.runtime.scriptrunner import get_script_run_ctx
    if get_script_run_ctx() is None:
        print("\n" + "="*70)
        print("⚠️  This is a Streamlit application.")
        print("="*70)
        print("\nPlease run it with the Streamlit CLI:")
        print(f"\n  streamlit run {os.path.abspath(__file__)}\n")
        print("="*70 + "\n")
        sys.exit(0)
except ImportError:
    pass

st.set_page_config(
    page_title="PubMed Genetic Variant Extractor",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 PubMed Genetic Variant & Clinical Data Extractor")

st.markdown("""
This tool searches PubMed publications for genetic variants and extracts clinical data using AI analysis.
Enter a gene name and optional variant to find relevant research papers.
""")

if 'results' not in st.session_state:
    st.session_state.results = None
if 'articles_analyzed' not in st.session_state:
    st.session_state.articles_analyzed = 0
if 'target_phenotypes' not in st.session_state:
    st.session_state.target_phenotypes = []

with st.sidebar:
    st.header("⚙️ Configuration")
    
    st.subheader("LLM Model Selection")
    model_choice = st.radio(
        "Choose AI Model:",
        ["gemini", "anthropic", "openai"],
        help="Gemini 2.5 Flash (Google - FREE tier), Claude 3 Haiku (Anthropic), or GPT-3.5 (OpenAI)"
    )
    
    api_key_available = False
    if model_choice == "gemini":
        if os.environ.get('GEMINI_API_KEY'):
            st.success("✓ Gemini API key detected")
            api_key_available = True
        else:
            st.warning("⚠️ GEMINI_API_KEY not set")
            st.info("Get a free API key from https://aistudio.google.com/apikey")
    elif model_choice == "anthropic":
        if os.environ.get('ANTHROPIC_API_KEY'):
            st.success("✓ Anthropic API key detected")
            api_key_available = True
        else:
            st.warning("⚠️ ANTHROPIC_API_KEY not set")
            st.info("Please add your API key in Secrets")
    else:
        if os.environ.get('OPENAI_API_KEY'):
            st.success("✓ OpenAI API key detected")
            api_key_available = True
        else:
            st.warning("⚠️ OPENAI_API_KEY not set")
            st.info("Please add your API key in Secrets")
    
    st.divider()
    
    st.subheader("Search Settings")
    max_results = st.slider(
        "Maximum articles to analyze:",
        min_value=5,
        max_value=200,
        value=50,
        step=5,
        help="Increase to analyze more articles (may be slower and use more API credits)"
    )

    phenotype_input = st.text_input(
        "Target phenotypes (comma separated):",
        value="sex, QT interval",
        help="List key phenotypes or clinical features for the AI to prioritize (e.g., sex, QT interval, arrhythmia)"
    )
    
    email = st.text_input(
        "Email (for PubMed API):",
        value="user@example.com",
        help="Required by NCBI for PubMed API access"
    )

col1, col2 = st.columns([2, 1])

with col1:
    gene_input = st.text_input(
        "Gene Name:",
        placeholder="e.g., KCNH2, BRCA1, TP53",
        help="Enter the gene symbol you want to search for"
    )

with col2:
    variant_input = st.text_input(
        "Variant (optional):",
        placeholder="e.g., p.Tyr54Asn",
        help="Specific variant notation (optional)"
    )

st.markdown("**Example searches:**")
st.markdown("- Gene: `KCNH2`, Variant: `p.Tyr54Asn` (QT prolongation)")
st.markdown("- Gene: `BRCA1`, Variant: `p.Arg1699Gln` (breast cancer)")
st.markdown("- Gene: `SCN5A` (leave variant empty for all variants)")

search_button = st.button("🔍 Search PubMed & Analyze", type="primary", disabled=not api_key_available)

if search_button:
    if not gene_input:
        st.error("Please enter a gene name")
    else:
        with st.spinner("Searching PubMed and analyzing articles with AI..."):
            try:
                pubmed = PubMedAPI(email=email)
                st.info(f"Searching PubMed for {gene_input}" + (f" {variant_input}" if variant_input else ""))
                
                articles = pubmed.search_and_fetch(
                    gene=gene_input,
                    variant=variant_input if variant_input else None,
                    max_results=max_results
                )
                
                if not articles:
                    st.warning("No articles found. Try different search terms.")
                else:
                    full_text_count = sum(1 for a in articles if a.get('content_type') == 'full_text_with_supplements')
                    abstract_count = len(articles) - full_text_count
                    
                    content_summary = []
                    if full_text_count > 0:
                        content_summary.append(f"{full_text_count} with full text & supplements")
                    if abstract_count > 0:
                        content_summary.append(f"{abstract_count} abstract-only")
                    
                    st.success(f"Found {len(articles)} articles ({', '.join(content_summary)}). Analyzing with AI...")
                    
                    progress_bar = st.progress(0)
                    status_text = st.empty()

                    analyzer = LLMAnalyzer(model_choice=model_choice)

                    target_phenotypes = [
                        p.strip() for p in phenotype_input.split(",") if p.strip()
                    ] if phenotype_input else []

                    results = []
                    for idx, article in enumerate(articles):
                        content_label = "full text" if article.get('content_type') == 'full_text_with_supplements' else "abstract"
                        status_text.text(f"Analyzing article {idx + 1} of {len(articles)} ({content_label})...")
                        progress_bar.progress((idx + 1) / len(articles))

                        analysis = analyzer.extract_variant_data(
                            article['full_text'],
                            gene_input,
                            variant_input if variant_input else None,
                            target_phenotypes=target_phenotypes
                        )
                        
                        result = {
                            "pmid": article['pmid'],
                            "title": article['title'],
                            "year": article['year'],
                            "authors": article['authors'],
                            "has_variant_data": analysis.get('has_variant_data', False),
                            "study_metadata": analysis.get('study_metadata', {}),
                            "individual_patients": analysis.get('individual_patients', []),
                            "aggregate_data": analysis.get('aggregate_data', []),
                            "confidence": analysis.get('confidence', 'low'),
                            "content_type": article.get('content_type', 'abstract'),
                            "pmc_id": article.get('pmc_id')
                        }

                        results.append(result)
                    
                    st.session_state.results = results
                    st.session_state.target_phenotypes = target_phenotypes
                    st.session_state.articles_analyzed = len(articles)
                    
                    progress_bar.empty()
                    status_text.empty()
                    st.success(f"✅ Analysis complete! Processed {len(articles)} articles.")
                    
            except Exception as e:
                st.error(f"Error: {str(e)}")

if st.session_state.results:
    st.divider()
    st.header("📊 Results")

    results = st.session_state.results
    target_phenotypes = st.session_state.get('target_phenotypes', [])

    if target_phenotypes:
        st.caption(
            "Focused phenotypes: " + ", ".join(target_phenotypes)
        )

    articles_with_data = sum(1 for r in results if r['has_variant_data'])
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Articles Analyzed", len(results))
    with col2:
        st.metric("Articles with Variant Data", articles_with_data)
    with col3:
        st.metric("Detection Rate", f"{(articles_with_data/len(results)*100):.1f}%")
    
    # Collect deduplication warnings
    cohort_names = {}
    for result in results:
        metadata = result.get('study_metadata', {})
        cohort_name = metadata.get('cohort_name')
        if cohort_name and cohort_name.lower() not in ['null', 'none', 'n/a', '']:
            if cohort_name not in cohort_names:
                cohort_names[cohort_name] = []
            cohort_names[cohort_name].append(f"PMID {result['pmid']} ({result['year']})")

    # Show deduplication warnings
    if cohort_names:
        duplicated_cohorts = {k: v for k, v in cohort_names.items() if len(v) > 1}
        if duplicated_cohorts:
            st.warning("⚠️ **Potential Participant Overlap Detected**")
            st.markdown("The following cohorts appear in multiple studies. Participants may be duplicated:")
            for cohort, pmids in duplicated_cohorts.items():
                st.markdown(f"- **{cohort}**: {', '.join(pmids)}")
            st.markdown("*Review these studies carefully to avoid double-counting participants.*")
            st.divider()

    # SECTION 1: Individual Patient Data
    st.subheader("👤 Individual Patient Data")
    st.caption("Case reports and individual patient narratives from publications")

    all_patients = []
    for result in results:
        if result['has_variant_data'] and result.get('individual_patients'):
            metadata = result.get('study_metadata', {})
            for patient in result['individual_patients']:
                all_patients.append({
                    'PMID': result['pmid'],
                    'Title': result['title'][:60] + '...' if len(result['title']) > 60 else result['title'],
                    'Year': result['year'],
                    'Study Type': metadata.get('study_type', 'N/A'),
                    'Location': metadata.get('location', 'N/A'),
                    'Cohort': metadata.get('cohort_name', 'N/A'),
                    'Gene': patient.get('gene', 'N/A'),
                    'Variant': patient.get('variant', 'N/A'),
                    'Genotype': patient.get('genotype', 'N/A'),
                    'Phenotypes': ', '.join(patient.get('phenotypes', [])) if patient.get('phenotypes') else 'N/A',
                    'Age': patient.get('age', 'N/A'),
                    'Sex': patient.get('sex', 'N/A'),
                    'Ethnicity': patient.get('ethnicity', 'N/A'),
                    'Treatment': patient.get('treatment', 'N/A'),
                    'Outcome': patient.get('outcome', 'N/A'),
                    'Confidence': result['confidence']
                })

    if all_patients:
        df_patients = pd.DataFrame(all_patients)
        st.dataframe(df_patients, use_container_width=True, hide_index=True)
        st.caption(f"📊 Total individual patients: {len(all_patients)}")
    else:
        st.info("No individual patient data found in analyzed articles.")

    st.divider()

    # SECTION 2: Aggregate Study Data
    st.subheader("📊 Aggregate Study Data (Case/Control Counts)")
    st.caption("Population-level statistics from cohort studies, GWAS, and meta-analyses")

    all_aggregate = []
    for result in results:
        if result['has_variant_data'] and result.get('aggregate_data'):
            metadata = result.get('study_metadata', {})
            for agg in result['aggregate_data']:
                all_aggregate.append({
                    'PMID': result['pmid'],
                    'Title': result['title'][:60] + '...' if len(result['title']) > 60 else result['title'],
                    'Year': result['year'],
                    'Study Type': metadata.get('study_type', 'N/A'),
                    'Location': metadata.get('location', 'N/A'),
                    'Cohort': metadata.get('cohort_name', 'N/A'),
                    'Gene': agg.get('gene', 'N/A'),
                    'Variant': agg.get('variant', 'N/A'),
                    'Cases w/ Variant': agg.get('cases_with_variant', 'N/A'),
                    'Total Cases': agg.get('total_cases', 'N/A'),
                    'Controls w/ Variant': agg.get('controls_with_variant', 'N/A'),
                    'Total Controls': agg.get('total_controls', 'N/A'),
                    'Ethnicity Breakdown': agg.get('ethnicity_breakdown', 'N/A'),
                    'AF Cases': agg.get('allele_frequency_cases', 'N/A'),
                    'AF Controls': agg.get('allele_frequency_controls', 'N/A'),
                    'Phenotype': agg.get('phenotype', 'N/A'),
                    'Notes': agg.get('notes', 'N/A'),
                    'Confidence': result['confidence']
                })

    if all_aggregate:
        df_aggregate = pd.DataFrame(all_aggregate)
        st.dataframe(df_aggregate, use_container_width=True, hide_index=True)
        st.caption(f"📊 Total aggregate variant entries: {len(all_aggregate)}")
    else:
        st.info("No aggregate study data found. Most articles may be case reports rather than population studies.")

    # Export section
    if all_patients or all_aggregate:
        st.divider()
        st.subheader("📥 Export Data")

        col1, col2 = st.columns(2)

        if all_patients:
            with col1:
                csv_patients = df_patients.to_csv(index=False)
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                filename_patients = f"individual_patients_{gene_input}_{timestamp}.csv"

                st.download_button(
                    label="Download Individual Patient Data",
                    data=csv_patients,
                    file_name=filename_patients,
                    mime="text/csv"
                )

        if all_aggregate:
            with col2:
                csv_aggregate = df_aggregate.to_csv(index=False)
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                filename_aggregate = f"aggregate_data_{gene_input}_{timestamp}.csv"

                st.download_button(
                    label="Download Aggregate Study Data",
                    data=csv_aggregate,
                    file_name=filename_aggregate,
                    mime="text/csv"
                )
    else:
        st.info("No variant data extracted from the analyzed articles. Try:")
        st.markdown("- Different gene or variant terms")
        st.markdown("- Increasing the number of articles")
        st.markdown("- Using more general search terms")
    
    with st.expander("📄 View All Article Details"):
        for result in results:
            content_badge = "📖 Full Text + Supplements" if result.get('content_type') == 'full_text_with_supplements' else "📄 Abstract Only"
            pmc_info = f" (PMC{result.get('pmc_id')})" if result.get('pmc_id') else ""

            st.markdown(f"**PMID {result['pmid']}{pmc_info}** - {result['title']}")
            st.markdown(f"*{result['authors']} ({result['year']})* | {content_badge}")
            st.markdown(f"Has variant data: {'✅ Yes' if result['has_variant_data'] else '❌ No'} | Confidence: {result['confidence']}")

            # Show study metadata
            metadata = result.get('study_metadata', {})
            if metadata and any(metadata.values()):
                st.markdown("**Study Metadata:**")
                st.json(metadata)

            # Show individual patients
            if result['has_variant_data'] and result.get('individual_patients'):
                st.markdown(f"**Individual Patients ({len(result['individual_patients'])}):**")
                for v in result['individual_patients']:
                    st.json(v)

            # Show aggregate data
            if result['has_variant_data'] and result.get('aggregate_data'):
                st.markdown(f"**Aggregate Data ({len(result['aggregate_data'])}):**")
                for a in result['aggregate_data']:
                    st.json(a)

            st.divider()

st.divider()
st.markdown("""
### How it works:
1. **Search PubMed**: Queries the PubMed database for articles related to your gene/variant
2. **Retrieve Full Text**: Automatically fetches full text + supplements from PubMed Central (PMC) when available, otherwise uses abstracts
3. **AI Analysis**: Uses AI to extract TWO types of data from articles:
   - **Individual Patient Data**: Case reports and patient narratives with demographics, phenotypes, genotype, treatment, and outcomes
   - **Aggregate Study Data**: Population-level statistics including case/control counts, ethnicity breakdowns, and allele frequencies
4. **Study Metadata**: Extracts study location, cohort names, and recruitment information to help identify potential participant overlap
5. **Deduplication Alerts**: Warns when the same cohort appears in multiple studies to avoid double-counting participants
6. **Export Results**: Download separate CSV files for individual patient data and aggregate study statistics

**New Features:**
- 📍 **Location & Ethnicity**: Tracks study geography and population ancestry
- 📊 **Case/Control Counts**: Extracts sample sizes from cohort studies and GWAS
- ⚠️ **Overlap Detection**: Identifies when studies use the same cohorts (e.g., UK Biobank, CHARGE Consortium)

**Note**: Full text with supplements is available for open-access articles in PMC. Other articles will be analyzed using abstracts only.
""")
