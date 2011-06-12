import json

def create_json(query,subject):
  "formats my file information into a json file needed for coannot"
  json_file = open("data/{0}_{1}/{0}_{1}.json".format(query,subject), "w")
  json_dict = {"genome_a": {"name":"{0}".format(query),"flat":"data/{0}_{1}/{0}.bed".format(query,subject),
  "fasta":"data/{0}_{1}/{0}.fasta".format(query,subject)},"genome_b": {"name":"{0}".format(subject),"flat":"data/{0}_{1}/{1}.bed".format(query,subject),
  "fasta":"data/{0}_{1}/{1}.fasta".format(query,subject)},"blast": {"W": 20,"a": 8,"e": 0.001},
  "defult ": {"blast_path":"/user/bin/blastall","out_dir":"data/{0}_{1}/".format(query,subject),"min_len":100, "reciprocal":True,"blast_log":True}}
  json_formatted = json.dumps(json_dict,indent=4)
  json_file.write(json_formatted)

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--query", dest="query", help="name of query org or ORGA")
    parser.add_option("--subject", dest="subject", help="name of subject org or ORGB")
    (options, _) = parser.parse_args()
    create_json('options.query', 'options.subject')
    #create_json(options.query, options.subject)

