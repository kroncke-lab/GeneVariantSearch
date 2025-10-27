import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pubmed_api import PubMedAPI
from llm_analyzer import LLMAnalyzer
import json

class VariantValidator:
    def __init__(self, api_key: str, model_choice: str):
        self.api_key = api_key
        self.model_choice = model_choice
        self.pubmed = PubMedAPI()
        self.analyzer = LLMAnalyzer(api_key, model_choice)
    
    def validate_variant(self, pmid: str, annotated_variant: str, annotated_phenotype: str) -> dict:
        """
        Validate a variant annotation against the source paper
        
        Returns:
            dict with keys: validation_status, validation_reason, paper_title, reported_data
        """
        try:
            articles = self.pubmed.search_pubmed(query=f"{pmid}[PMID]", max_results=1)
            
            if not articles:
                return {
                    'pmid': pmid,
                    'validation_status': 'UNSURE',
                    'validation_reason': 'Paper not found in PubMed',
                    'paper_title': '',
                    'reported_data': ''
                }
            
            article = articles[0]
            paper_title = article.get('title', '')
            
            text_content = f"Title: {paper_title}\n\nAbstract: {article.get('abstract', '')}"
            
            pmc_id = article.get('pmc_id')
            if pmc_id:
                fulltext_data = self.pubmed.fetch_pmc_fulltext(pmc_id, gene=annotated_variant.split(':')[0] if ':' in annotated_variant else '')
                if fulltext_data.get('full_text'):
                    text_content += f"\n\nFull Text:\n{fulltext_data['full_text'][:20000]}"
                if fulltext_data.get('supplements'):
                    text_content += f"\n\nSupplemental Data:\n{fulltext_data['supplements'][:30000]}"
            
            validation_prompt = f"""You are validating genetic variant annotations against the source publication.

**Annotated Information from CSV:**
- Variant: {annotated_variant}
- Phenotype/Diagnosis: {annotated_phenotype}

**Source Paper Content:**
{text_content[:50000]}

**Task:**
Compare the annotated information against what is actually reported in the paper. Determine:
1. Is the variant correctly identified? (check variant name, HGVS notation, gene)
2. Is the phenotype/diagnosis correctly annotated? (check reported clinical features, disease, outcomes)
3. Overall assessment: CORRECT, INCORRECT, or UNSURE

**Output JSON format:**
{{
    "validation_status": "CORRECT" | "INCORRECT" | "UNSURE",
    "variant_match": "exact" | "partial" | "mismatch" | "not_found",
    "phenotype_match": "exact" | "partial" | "mismatch" | "not_found",
    "reported_variant": "what the paper actually says about the variant",
    "reported_phenotype": "what the paper actually says about phenotype/diagnosis",
    "validation_reason": "detailed explanation of your assessment"
}}

**Validation Guidelines:**
- CORRECT: Both variant and phenotype match the paper (exact or semantically equivalent)
- INCORRECT: Clear mismatch between annotation and paper (wrong variant, wrong phenotype, or contradictory info)
- UNSURE: Paper doesn't provide enough information, ambiguous reporting, or unclear match
- Be strict but fair - minor notation differences (p. vs protein name) are OK if semantically identical
- Consider that abbreviations and full terms may be used interchangeably
"""

            response = self.analyzer._call_llm(validation_prompt)
            
            try:
                result_data = json.loads(response)
            except json.JSONDecodeError:
                json_match = response.find('{')
                if json_match != -1:
                    json_end = response.rfind('}') + 1
                    result_data = json.loads(response[json_match:json_end])
                else:
                    result_data = {
                        'validation_status': 'UNSURE',
                        'validation_reason': 'AI response parsing error',
                        'reported_variant': '',
                        'reported_phenotype': ''
                    }
            
            return {
                'pmid': pmid,
                'validation_status': result_data.get('validation_status', 'UNSURE'),
                'validation_reason': result_data.get('validation_reason', 'No reason provided'),
                'paper_title': paper_title,
                'reported_data': f"Variant: {result_data.get('reported_variant', 'N/A')} | Phenotype: {result_data.get('reported_phenotype', 'N/A')}"
            }
            
        except Exception as e:
            return {
                'pmid': pmid,
                'validation_status': 'ERROR',
                'validation_reason': f'Validation error: {str(e)}',
                'paper_title': '',
                'reported_data': ''
            }
