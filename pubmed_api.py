import os
import time
from Bio import Entrez
import xml.etree.ElementTree as ET
from typing import List, Dict, Optional

Entrez.email = os.environ.get("PUBMED_EMAIL", "user@example.com")

class PubMedAPI:
    def __init__(self, email: str = None):
        if email:
            Entrez.email = email
        self.max_results = 20
    
    def search_publications(self, gene: str, variant: Optional[str] = None, max_results: int = 20) -> List[str]:
        query = gene
        if variant:
            query = f"{gene} AND {variant}"
        
        query += " AND (genetic variant OR mutation OR polymorphism OR genotype)"
        
        try:
            search_handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=max_results,
                sort="relevance"
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            id_list = search_results.get("IdList", [])
            return id_list
        except Exception as e:
            raise Exception(f"PubMed search failed: {str(e)}")
    
    def fetch_article_details(self, pmid_list: List[str]) -> List[Dict]:
        if not pmid_list:
            return []
        
        articles = []
        
        try:
            fetch_handle = Entrez.efetch(
                db="pubmed",
                id=pmid_list,
                rettype="xml",
                retmode="xml"
            )
            
            xml_data = fetch_handle.read()
            fetch_handle.close()
            
            root = ET.fromstring(xml_data)
            
            for article in root.findall(".//PubmedArticle"):
                try:
                    pmid_elem = article.find(".//PMID")
                    pmid = pmid_elem.text if pmid_elem is not None else "Unknown"
                    
                    title_elem = article.find(".//ArticleTitle")
                    title = title_elem.text if title_elem is not None else "No title"
                    
                    abstract_texts = article.findall(".//AbstractText")
                    abstract = " ".join([
                        elem.text for elem in abstract_texts if elem.text
                    ]) if abstract_texts else "No abstract available"
                    
                    year_elem = article.find(".//PubDate/Year")
                    year = year_elem.text if year_elem is not None else "Unknown"
                    
                    authors = []
                    for author in article.findall(".//Author"):
                        lastname = author.find("LastName")
                        forename = author.find("ForeName")
                        if lastname is not None and forename is not None:
                            authors.append(f"{forename.text} {lastname.text}")
                    
                    author_str = ", ".join(authors[:3]) if authors else "Unknown authors"
                    if len(authors) > 3:
                        author_str += " et al."
                    
                    articles.append({
                        "pmid": pmid,
                        "title": title,
                        "abstract": abstract,
                        "year": year,
                        "authors": author_str,
                        "full_text": f"{title}\n\n{abstract}"
                    })
                    
                except Exception as e:
                    continue
            
            return articles
            
        except Exception as e:
            raise Exception(f"Failed to fetch article details: {str(e)}")
    
    def search_and_fetch(self, gene: str, variant: Optional[str] = None, max_results: int = 20) -> List[Dict]:
        pmid_list = self.search_publications(gene, variant, max_results)
        
        if not pmid_list:
            return []
        
        time.sleep(0.34)
        
        articles = self.fetch_article_details(pmid_list)
        return articles
