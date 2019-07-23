import argparse
import csv
import os
import re
import sys
from textwrap import wrap


#For each entry in GenBank file retrieves nucleotide sequence, collection date and country 
#Generates file with sequences in fasta format, 
#sequences' names are in the following format: GenbankAccessionNumber_country_collectionYear
#saves the file in the directory of input_file

#input_file - file with sequences in GenBank format
#min length, max length - minimal and maximal lenghts of sequence in origin field
#type - type of alignment, short of full fragment of gene
#output - name of output file
def parse_gb(input_file, min_length, max_length):

    #File with countries/cities/other features and their abbreviations
    #Should be located at the same directory as script
    COUNTRY_MAP_FILE = os.path.join(sys.path[0],"country_map.csv")
    FEATURE_MAP_FILE = os.path.join(sys.path[0],"feature_map.csv") #not developed yet
    CITY_MAP_FILE = os.path.join(sys.path[0],"city_map.csv")

    #name of output file

    OUTPUT_FILE = '.'.join(input_file.split('.')[:-1]) + '.fasta'


    #length tresholds for sequences extracted from origin field

    MIN_ORIGIN_SIZE = min_length
    MAX_ORIGIN_SIZE = max_length

    #fields in GB entry which will be extracted
    FIELD_NAMES = ["strain", "isolate", "country", "collection_date"]
    
    #country_map - file with abbreviations of countries
    country_map = read_csv(COUNTRY_MAP_FILE)
    
    #feature_map = read_csv(FEATURE_MAP_FILE, strip_it=False)
    #city_map = read_csv(CITY_MAP_FILE)
    
    #RegExp for parsing  INPUT_FILE
    accession = re.compile(r"^ACCESSION\s+([a-z_A-Z0-9]+)") #accession number
    features = re.compile(r"^\s+\/([a-z_]+)=\"([^\"]+)\"$") #qualifiers from FEATURE source field

    #origin = re.compile(r"\s+\d+\s+([a-z ]+)$")
    origin = re.compile(r"\s+\d+\s+([atgcnrykmswbdhvu\s]+)$") #origin field (contains nt sequence)
    end_line = re.compile(r"^//\s+$") #end line of entry
    clean_from = re.compile(r"[:;.,/\s]")
    clean_to = "-"

    #
    def check_year(stri):
        year0 =  re.compile(r"[0-9]{4}")
        #_1994_ 
        year1 = re.compile(r"[_\"\-\s,;/]+[0-9]{4}[_\"\-\s,;]+")
        m = year0.search(stri)
        if m:
            return year0.search(m.group()).group()


    #Parsing INPUT_FILE

    tests = [] # stores entries
    tests_nodate = [] # list for entries with no collection date

    test_accession = "" # accession number of entry
    test_features = {}  #dictionary for qualifiers from FEATURE source field
    test_origin = "" #ORIGIN field
    
    if not os.path.exists(input_file):
        print('{filename}: No such file or directory'.format(filename=input_file))
        return


    with open(input_file, "r") as in_f:
        for line in in_f:

            m = accession.match(line) #finds ACCESSION field using RegExp
            if (m):
                test_accession = m.group(1)
                
            #check whether /source has started or has ended
            if re.match(r"^\s{5}[a-zA-z_]+[\s0-9<>\.]+", line) or re.match(r"^[a-zA-z]+", line):
                s = 0 # we are outside source field
            if re.match(r"^\s+source[0-9\.\s]+", line):
                s = 1 #inside source field
    
            if s == 1:
                m = features.match(line) #finds all features in field 'FEATURES source' using RegExp
                if (m):
                    test_features[m.group(1)] = m.group(2)

            m = origin.match(line) #ORIGIN field  (contains na sequence)
            if(m):
                test_origin += re.sub("\s+","",m.group(1))


            if(end_line.search(line)): #the end of entry

                if ("collection_date" or "collected_by") in test_features: 
                    test_features["collection_date"] = check_year(test_features["collection_date"]) #reformates collection date

                if "country" in test_features:
                    test_features["country"] = map_feature(test_features["country"], country_map) #replaces country by its abbreviation


                #for k in ["strain", "organism"]:
                    #if k in test_features: 
                        #test_features[k] = map_feature(test_features[k], feature_map)
            
            
                if ("collection_date" or "collected_by") not in test_features: #adds entries without collection date to test_nodate list
                    tests_nodate.append({
                        "accession" : test_accession, 
                        "features" : test_features,
                        "origin" : test_origin
                        })
                
                else:
                    tests.append({
                        "accession" : test_accession, 
                        "features" : test_features,
                        "origin" : test_origin
                        })

                #resets variables
                test_accession = "" 
                test_features = {}
                test_origin = ""
            
    in_f.close()

    #adds entries without collection date to the end of tests list
    for ent in tests_nodate:
        tests.append(ent)

    # Writing sequences in fasta-format to OUTPUT_FILE 
    out = open(OUTPUT_FILE, "w+")
    for test in tests:
        aaccession, f, origin = test["accession"], test["features"], test["origin"]
        #print(aaccession, f, origin)
        #print(aaccession)
        if not (MIN_ORIGIN_SIZE < len(origin) < MAX_ORIGIN_SIZE):
            print(aaccession, len(origin))
            continue

        #new name of entry = GB accession number + country + year
        output_name = ">%s" % aaccession
        for field_name in FIELD_NAMES:
            if field_name not in f:
                continue
            output_name += "_%s" % re.sub(clean_from, clean_to,f[field_name])
        output_name += "\n"
        out.write(output_name)
        
        #writes sequence from ORIGIN field
        out.write('\n'.join(wrap(origin,60)) + '\n')
        
    out.close()
    return OUTPUT_FILE

#custom csv parser :/
def read_csv(file_name, strip_it=True):
    if not os.path.exists(file_name):
        return {}

    def strip(value):
        return strip_it and value.strip() or value

    with open(file_name) as csvfile:
        reader = csv.DictReader(csvfile,
                                delimiter=",", 
                                fieldnames=["base", "new"])
        result = {}
        for row in reader:
            result[strip(row["base"].lower())] = strip(row["new"])
        csvfile.close()
        return result

#feature_map - dictionary, e.g. feature_map['Italia']='ITA'
#feature - possible key from feature_map
#return value if key is in feature_map
def map_feature(feature, feature_map):
    for k, v in feature_map.items():
        if feature.lower().startswith(k):
            return v
    return feature

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required=True)
    parser.add_argument("-min", "--min_length", type=int,
                        help="Minimal length of sequence.\
                        Sequences shorter than min length will not be included in the final dataset",\
                        required=True)
    parser.add_argument("-max", "--max_length", type=int,
                        help="Maximal length of sequence. \
                        Sequences longer than max length will not be included in the final dataset",\
                        required=True)
    args = parser.parse_args()

    if not len(sys.argv) == 7:
        print("Please, use \"python parser_gb.py --help\"")
    else:
        parse_gb(args.input_file, args.min_length, args.max_length)
