use std::error::Error;
use crate::{utils::Text, alphabets::dna};
use reqwest;
use std::collections::HashMap;
use serde_json::{Value};

static URL_NCBI: &'static str = "https://bast.ncbi.nlm.nih.gov/Blast.cgi";

pub fn create_fasta_from_adn(dna_str : Text) -> String {
    let mut dna_fasta = String::from(">Search dna from rust-bio \n");
    let mut count = 0;
    for i in dna_str{
        let nc = i as char;
        if count < 80{
            dna_fasta.push(nc);
        }
        else {
            count = 0;
            dna_fasta.push('\n');
            dna_fasta.push(nc);
        }
    }

    dna_fasta.to_owned()
}


pub async fn send_http_req(dna_fasta : String) -> Result<String,Box<dyn Error>> {
    let client_ncbi = reqwest::Client::builder().build()?;
    let params = [
        (String::from("QUERY"),dna_fasta),
        (String::from("DATABASE"),String::from("nt")),
        (String::from("PROGRAM"),String::from("blastn")),
        (String::from("CMD"),String::from("Put"))
        ];

    let resp = client_ncbi.put(URL_NCBI)
    .form(&params).send().await?;
    let resp_text = resp.text().await?;
    println!("{:#?}", resp_text);

    let index_rid = resp_text.find("RID = ").unwrap();
    let start_rid = index_rid+6;
    let iend_rid= resp_text.find("RTOE = ").unwrap();
    let end_rid = iend_rid-1;
    let rid = &resp_text[start_rid..end_rid];


    let params_result = [
        ("CMD","Get"),
        ("RID",rid)
    ];
    let mut status = String::from("WAITI");
    let mut search_resp = client_ncbi.get(URL_NCBI)
    .form(&params_result).send().await?;

    while status != "READY"{
        search_resp = client_ncbi.get(URL_NCBI)
        .form(&params_result).send().await?;
        let resp_search_text = search_resp.text().await?;
        let index_status = resp_search_text.find("Status=").unwrap() +7;
        let max_index = index_status + 5;
        status.clear();
        status = (&resp_search_text[index_status..max_index]).to_string();
        
    }

    let json_params_result = [
        ("CMD","Get"),
        ("RID",rid),
        ("FORMAT_TYPE","JSON2_S")
    ];
    search_resp = client_ncbi.get(URL_NCBI)
        .form(&json_params_result).send().await?;
    let json_text = search_resp.text().await?;

    let parsed = read_json(&json_text);
    let list_res = parsed["hits"][0].clone();
    let blast_result = list_res.as_str().unwrap().to_string();


    Ok(blast_result)
}

fn read_json(raw_json : &str) -> Value {
    let parsed: Value = serde_json::from_str(raw_json).unwrap();
    parsed
}

pub fn blast_req(){
    
}

mod tests {
    use crate::blast::blast_req::send_http_req;
    //use super::*;

    #[test]
    fn print_html_response() {
        let fasta_example = String::from(">Ejemplo de algo
        ATGGTAAGTTATTCTTCTAGTTGATATATTCTTACTCTTTCTTTCTTACGTAACAACTGATCAAGTGTTA
        ATATATTGTATTTAATCAGCAAGGCCAGTGGATAGCCGCAAGAGACCTTTCTATTACATGGGTGGACAAT
        CCTCAGTACTGGACATGGAAAACTGTTGATCCTAAGTCTGTCTCTCTCTAACTATCTCTCTCTAAATTTC
        ATTAAAAAGAGAAATATATTTTAATCCGATCTTAACATTTTGTACTACAGTATTGAAGTGGCAGAGCTTC
        TTAGGGTAGCTTGGCTTGACATTTATGGAAAGATCGAGACAAAAAATCTTATTCGAAAGACTAGTTATGC
        TGTATATTTAGTGTTCAAGTTAACAGATAACCCTCGTGAACTTGAACGAGCCACAGCGTCGCTAAGATTT
        GTGAACGAAGTGGCGGAGGGCGCTGGCATTGAGGGTACCACTGTTTTCATCTCGAAGAAAAAGGAATTAC
        CAGGAGAACTTGGCCGGTTCCCACATCTCCGAAGTGATGGCTGGTTAGAAACCAAGCTTGGTGAGTTTTT
        CAACAACTTAGGAGAGGATGGTGAAGTCGAAATGAGGTTGATGGAAATCAATGACAAAACTTGGAAATCT
        GGCATCATTGTTAAGGGCTTCGACATTCGTCCAAACTAA");
        let _result = send_http_req(fasta_example);
        
        
    }
}
/*
1 Func e volver adn a fasta DONE
2 generar el request
3 blast req principal publico
4 aplicarlo a blast p tamb
 */