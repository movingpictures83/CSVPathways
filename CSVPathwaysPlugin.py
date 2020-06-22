import sys
import pythoncyc.config as config
import pythoncyc
#this creates PGDB object associated with self.meta(MetaCyc)

from pythoncyc.PToolsFrame import PFrame


def findNonDisjointSets(sets):
   for i in range(len(sets)):
      for j in range(i+1, len(sets)):
         if (not sets[i].isdisjoint(sets[j])):
            return [i, j]
   return -1

class CSVPathwaysPlugin:
   def input(self, inputfile):
      params = open(inputfile, 'r')
      self.parameters = dict()
      for line in params:
         contents = line.strip().split('\t')
         self.parameters[contents[0]] = contents[1]
      config.set_host_name('castalia.cs.fiu.edu')
      self.meta = pythoncyc.select_organism('meta')
      abundancefile = open(self.parameters['csvfile'], 'r')
      mappingfile = open(self.parameters['mapping'], 'r')

      junk = mappingfile.readline()
      mapping = dict()
      self.reversemapping = dict()
      for line in mappingfile:
         myline = line.strip()
         contents = myline.split('\t')
         mapping["X"+contents[0]] = contents[1]
         if (contents[1] != "NOTFOUND"):
            self.reversemapping[contents[1]] = "X"+contents[0]
 

      firstline = abundancefile.readline().strip()
      entries = firstline.split(',')
      entries.remove(entries[0])
      self.microbes = dict()
      self.metabolites = []
      #print(entries)
      for entry in entries:
         # Microbe
         entry2 = entry[1:len(entry)-1]
         if entry2[len(entry2)-3] == '.' or entry2 == "Unassigned":
            # Not classified as lowest level
            if entry2.find('.') != len(entry2)-3:
               microbe = entry2[2:len(entry2)-3]
            else:
               microbe = entry2[0:len(entry2)-3]
            self.microbes[microbe] = entry
         # Metabolite
         else:
            metabolite = mapping[entry2]
            if (metabolite != "NOTFOUND"):
               self.metabolites.append(metabolite)

      #print(self.microbes)
      #print(self.metabolites)
      #raw_input()

   def run(self):
      
      self.sets = []

      self.pathwayset = set()

      for metabolite in self.metabolites:
         print("Computing pathways for: "+str(metabolite))
         try:
           for pathway in self.meta.pathways_of_compound(metabolite):
      
      
            x = set()
      
            flag = False
            if ('species' in PFrame(pathway, self.meta, getFrameData=True).__dict__):
             for species in PFrame(pathway, self.meta, getFrameData=True).__dict__['species']:
               if ("human" in PFrame(species, self.meta, getFrameData=True).__dict__['names']):
                  flag = True
               for name in PFrame(species, self.meta, getFrameData=True).__dict__['names']:
                  for microbe in self.microbes.keys():
                     if name.find(microbe) != -1:
                        x.add(self.microbes[microbe])
                        flag = True
            # Keep this path if it (1) occurs in a human and/or (2) contains microbe(s) from our sample
            if (flag):
               self.pathwayset.add(pathway)
               for compound in self.meta.compounds_of_pathway(PFrame(pathway, self.meta, getFrameData=True).__dict__['frameid']):
                   if (compound[1:len(compound)-1] in self.reversemapping):
                     x.add(self.reversemapping[compound[1:len(compound)-1]])
      
            #if (len(x) > 0):
            #   flag2 = False
            #   for entry in self.sets:
            #      if not entry.isdisjoint(x):
            #         self.sets.remove(entry)
            #         self.sets.append(entry.union(x))
            #         flag2 = True
            #         break
            #   if not flag2:
            self.sets.append(x)
         except:
            print("Warning: Mapping failed for compound "+str(metabolite)+".")
            print("Error was: "+str(sys.exc_info()))
            continue

   def output(self, outputfile):
            
      txtfile = open(outputfile+".txt", 'w')
      # Output txt file
      for pway in self.pathwayset:
         y = PFrame(pway, self.meta, getFrameData=True).__dict__ 
         #txtfile.write("*********************************************************************************************\n")
         txtfile.write("PATHWAY: "+y['common_name'].replace("<i>","").replace("</i>","").replace("<b>","").replace("</b>",""))
         #txtfile.write("\n")
         #try:
         #   txtfile.write("COMMENT: "+y['comment'][0].replace("<i>","").replace("</i>","").replace("<b>","").replace("</b>",""))
         #except:
         #   txtfile.write("COMMENT: See database.")
         #txtfile.write("\n")
         txtfile.write(" INVOLVES: ")
         for compound in self.meta.compounds_of_pathway(y['frameid']):
               if (compound[1:len(compound)-1] in self.reversemapping):
                  #txtfile.write(compound+"\t")
                  #print compound[1:len(compound)-1], self.reversemapping[compound[1:len(compound)-1]]
                  txtfile.write(self.reversemapping[compound[1:len(compound)-1]]+"\t")
         if ('species' in PFrame(pway, self.meta, getFrameData=True).__dict__):
          #txtfile.write("\nOBSERVED IN: ")
          mfound = set()
          for species in PFrame(pway, self.meta, getFrameData=True).__dict__['species']:
                if ("human" in PFrame(species, self.meta, getFrameData=True).__dict__['names']):
                   mfound.add("human")
                   #txtfile.write(str(count)+". Human\n")
                else:
                   for name in PFrame(species, self.meta, getFrameData=True).__dict__['names']:
                      for microbe in self.microbes.keys():
                         if name.find(microbe) != -1:
                            mfound.add(microbe)
                            #break
          count = 1 
          for microbe in mfound:
               txtfile.write(microbe+"\t")
               #txtfile.write(str(count)+". "+microbe+"\n")
               count += 1
         txtfile.write("\n")
         #txtfile.write("*********************************************************************************************\n")
      
      # Merge self.sets
      ms = 0
      while (ms != -1):
         ms = findNonDisjointSets(self.sets)
         if (ms != -1):
            set1 = self.sets[ms[0]]
            set2 = self.sets[ms[1]]
            self.sets.remove(set1)
            self.sets.remove(set2)
            self.sets.append(set1.union(set2))
            
      
      #print(self.sets)
      count = 1
      noafile = open(outputfile+".noa", 'w')
      noafile.write('Name\tPathway\n')
      for myset in self.sets:
         for element in myset:
            noafile.write(element+"\t"+str(count)+"\n")
         count += 1
