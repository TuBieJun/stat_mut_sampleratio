extern crate clap;

use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::io::BufWriter;
use std::collections::HashMap;
use std::collections::HashSet;
use clap::{Arg, App, SubCommand};

fn parse_vcf(vcf_l:&str, panel:&str, vcf_t:&str,
            min_support:u32, min_fre:f64, min_depth:u32, 
            D_germline:&mut HashMap<String, HashMap<String, u32>>) -> u32 {
    
    let mut sample_num:u32 = 0;
    let mut pos_set:HashSet<String> = HashSet::new();

    let mut f = File::open(vcf_l)
                .expect("can't open vcf list file");
    let reader = BufReader::new(f); 

    for line in reader.lines() {
        let vcf_file = line.unwrap();
        let vcf_f = File::open(vcf_file);
        let vcf_f = match vcf_f {
            Ok(file) => file,
            Err(error) => {
                println!("There was a problem opening the file: {:?}", error);
                continue; 
            },
        };
        sample_num += 1;

        let reader_vcf = BufReader::new(vcf_f);
        for vcf_line in reader_vcf.lines() {
            let vcf_info = vcf_line.unwrap();
            if !vcf_info.starts_with("#") {
                
                let info_l:Vec<&str> = vcf_info.trim().split('\t').collect();
                let mut dp_info:Vec<&str>;
                let mut depth:u32 = 0;
                let mut support:u32 = 0;
                let mut fre:f64 = 0f64;
                if vcf_t == "gatk" {

                    dp_info = info_l[9].split(':').collect();
                    depth = dp_info[3].parse().unwrap();
                    let support_l:Vec<&str> = dp_info[1].split(',').collect();
                    support = support_l[1].parse().unwrap();
                    if support == 0 {
                        fre = 0f64;
                    } else {
                        fre = (support as f64)/(depth as f64);
                    }

                } else if vcf_t == "mutect" {

                    let mut n_index:usize;
                    if info_l[9].contains("0/1") {
                        n_index = 10;
                    }
                    else {
                        n_index = 9;
                    }
                    dp_info = info_l[n_index].split(':').collect();
                    depth = dp_info[3].parse().unwrap();
                    fre = dp_info[4].parse().unwrap();
                    let temp:Vec<&str> = dp_info[1].split(",").collect();
                    support = temp[0].parse().unwrap();

                }

                if support >= min_support && fre >= min_fre && depth >= min_depth {
                    let panel_s = String::from(panel);
                    let k = format!("{}\t{}\t{}\t{}", info_l[0], info_l[1], 
                                    info_l[3], info_l[4]);
                    let k1 = k.clone();
                    pos_set.insert(k1);
                    if D_germline.contains_key(&k) {
                        *D_germline.get_mut(&k).unwrap()
                            .entry(panel_s).or_insert(0) += 1;
                    }
                    else {
                        let mut panel_d = HashMap::new();
                        panel_d.insert(panel_s, 1);
                        D_germline.insert(k, panel_d);
                    }
                }
            }

        }
    }
    sample_num
}

fn out_info(out_file:&str, panel_list:&Vec<String>, D_germline:&HashMap<String, HashMap<String, u32>>, sample_num_list:&Vec<u32>) {
    let f = File::create(out_file).expect("can't create out file");
    let mut writer = BufWriter::new(f);
    let header = format!("chrom\tpos\tref\talt\t{}\n", panel_list.join("\t"));
    writer.write(header.as_bytes());
    //println!("{:?}", D_germline);
    for (k, v) in D_germline {
        let mut temp:Vec<String> = Vec::new();
        for (i,panel) in panel_list.iter().enumerate() {
           if v.contains_key(panel) {
               temp.push(format!("{}", (*v.get(panel).unwrap() as f64 / sample_num_list[i] as f64)));
           } else {
               temp.push(String::from("0"));
           }
        }
        let info = format!("{}\t{}\n", k, temp.join("\t"));
        writer.write(info.as_bytes());
    }

}



fn main() {
    let matches = App::new("stat the specific mut pos ratio")
                .version("1.0")
                .author("Teng Li <707704459@qq.com>")
                .arg(Arg::with_name("panel509").long("panel509").value_name("panel509")
                     .help("the input panel509 vcf list file"))
                .arg(Arg::with_name("exome").long("exome").value_name("exome")
                     .help("the input exome vcf list file"))
                .arg(Arg::with_name("wgs").long("wgs").value_name("wgs")
                     .help("the input wgs vcf list file"))
                .arg(Arg::with_name("out_file").short("o").long("out_file")
                     .value_name("out_file")
                     .help("the out file"))
                .arg(Arg::with_name("min_support").short("s")
                     .value_name("MIN_SUPPORT")
                     .default_value("4")
                     .help("the min specific mut pos support num, default is 4"))
                .arg(Arg::with_name("min_fre").short("f")
                     .value_name("MIN_FRE")
                     .default_value("0.2")
                     .help("the min specific mut pos support fre, default is 0.2"))
                .arg(Arg::with_name("min_depth").short("d")
                     .value_name("MIN_DEPTH")
                     .default_value("15")
                     .help("the min specific mut pos support depth, default is 15"))
                .get_matches();

    let min_support:u32 = matches.value_of("min_support").unwrap().trim().parse()
                        .expect("the -s arg must be a int number");
    let min_fre:f64 = matches.value_of("min_fre").unwrap().trim().parse()
                        .expect("the -f arg must be a float number");
    let min_depth:u32 = matches.value_of("min_depth").unwrap().trim().parse() 
                        .expect("the -d arg must be a int  number");
    let out_file = matches.value_of("out_file").unwrap();
    let mut panel_list:Vec<String> = Vec::new();
    let mut sample_num_list:Vec<u32> = Vec::new();
    let mut D_germline:HashMap<String, HashMap<String, u32>> = HashMap::new();
        
    
    if matches.is_present("panel509") {
        let panel509_list = matches.value_of("panel509").unwrap();
        let sample_num = parse_vcf(&panel509_list, "panel509", "mutect", 
                  min_support, min_fre, min_depth, &mut D_germline);
        panel_list.push(String::from("panel509"));
        sample_num_list.push(sample_num)
    }

    if matches.is_present("exome") {
        let exome_list = matches.value_of("exome").unwrap();
        let sample_num = parse_vcf(&exome_list, "exome", "mutect", 
                  min_support, min_fre, min_depth, &mut D_germline);
        panel_list.push(String::from("exome"));
        sample_num_list.push(sample_num)
    }

    if matches.is_present("wgs") {
        let wgs_list = matches.value_of("wgs").unwrap();
        let sample_num = parse_vcf(&wgs_list, "wgs", "gatk", 
                  min_support, min_fre, min_depth, &mut D_germline);
        panel_list.push(String::from("wgs"));
        sample_num_list.push(sample_num)
    }

    out_info(&out_file, &panel_list, &D_germline, &sample_num_list);
}
