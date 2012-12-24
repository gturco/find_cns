import json

def create_json(query,subject,blast):
  "formats my file information into a json file needed for coannot"
  json_file = open("data/{0}_{1}/{0}_{1}.json".format(query,subject), "w")
  json_dict = {"genome_a": {"name":"{0}".format(query),"flat":"data/{0}_{1}/{0}.bed".format(query,subject),
  "fasta":"data/{0}_{1}/{0}.fasta".format(query,subject)},"genome_b": {"name":"{0}".format(subject),"flat":"data/{0}_{1}/{1}.bed".format(query,subject),
  "fasta":"data/{0}_{1}/{1}.fasta".format(query,subject)},"blast": {"W": 20,"a": 8,"e": 0.001},
  "default": {"blast_path":"{0}/blastall".format(blast),"out_dir":"data/{0}_{1}/".format(query,subject),"min_len":100, "reciprocal":True,"blast_log":True}}
  json_formatted = json.dumps(json_dict,indent=4)
  json_file.write(json_formatted)

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--query", dest="query", help="name of query org or ORGA")
    parser.add_option("--subject", dest="subject", help="name of subject org or ORGB")
    parser.add_option("--blast_path",dest="blast", help="location of most up todate blast dir")    
    (options, _) = parser.parse_args()
    
    #try: open("data/{0}_{1}/{0}_{1}.json".format(options.query,options.subject), "w"):
    #except NameError:
    #    print "cannot  find dir {0}_{1} make sure to store data in {0}_{1}".format(options.query,options.subject)    
    create_json(options.query, options.subject,options.blast)

