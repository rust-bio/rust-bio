use std::error::Error;
use crate::utils::Text;
use reqwest;
use serde_json::{Value};
use std::fmt;

#[derive(Debug)]
struct MyError(String);

impl fmt::Display for MyError{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "There is an error: \n {}", self.0)
    }
}

impl Error for MyError {}

static URL_NCBI: &'static str = "https://bast.ncbi.nlm.nih.gov/Blast.cgi";

pub fn create_fasta_from_dna(dna_str : Text) -> String {
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
    let client_ncbi = reqwest::Client::builder().build()
        .or(Err(Box::new(MyError("Client not connected".into()))))?;
    let params = [
        (String::from("QUERY"),dna_fasta),
        (String::from("DATABASE"),String::from("nt")),
        (String::from("PROGRAM"),String::from("blastn")),
        (String::from("CMD"),String::from("Put"))
        ];

    let resp = client_ncbi.put(URL_NCBI)
        .form(&params).send().await
        .or(Err(Box::new(MyError("No response was recieved from client".into()))))?;
    let resp_text = resp.text().await
        .or(Err(Box::new(MyError("Response was not able to convert to text".into()))))?;
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
        .form(&params_result).send().await
        .or(Err(Box::new(MyError("No response was recieved from client with the ID".into()))))?;

    while status != "READY"{
        search_resp = client_ncbi.get(URL_NCBI)
            .form(&params_result).send().await
            .or(Err(Box::new(MyError("No response was recieved from client with the ID while waiting".into()))))?;
        let resp_search_text = search_resp.text().await
            .or(Err(Box::new(MyError("Response was not able to convert to text while waiting".into()))))?;
        let index_status = resp_search_text.find("Status=").unwrap() +7;
        let max_index = index_status + 5;
        status.clear();
        status = (&resp_search_text[index_status..max_index]).to_string();
        if status == "UNKWO"{
            return Err(Box::new(MyError("The ID saved was not found on the server".into())));
        }
        else{continue;}
        
    }

    let json_params_result = [
        ("CMD","Get"),
        ("RID",rid),
        ("FORMAT_TYPE","JSON2_S")
    ];
    search_resp = client_ncbi.get(URL_NCBI)
        .form(&json_params_result).send().await
        .or(Err(Box::new(MyError("Error on the response requested as JSON format".into()))))?;
    let json_text = search_resp.text().await
        .or(Err(Box::new(MyError("Error converting response JSON response to text".into()))))?;

    let parsed = read_json(&json_text);
    let list_res = parsed["hits"][0].clone();
    let blast_result = list_res.as_str().unwrap().to_string();

    if blast_result.is_empty() {
        Err(Box::new(MyError("No match found".into())))?
    }else{
    Ok(blast_result)}
}

fn read_json(raw_json : &str) -> Value {
    let parsed: Value = serde_json::from_str(raw_json).unwrap();
    parsed
}

pub async fn blastn_req(dna_str : Text) -> String {
    let dna_fasta = create_fasta_from_dna(dna_str);
    let results_blastn = send_http_req(dna_fasta).await;
    match results_blastn {
        Ok(n) => {return n},
        Err(e) => {
            panic!("Error found {}", e)
        },
    }
    
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