﻿using Proteomics;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class ProteinParsimonyEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private readonly bool treatModPeptidesAsDifferentPeptides;
        //dictionary with key being object contianing mass of the compact peptide and value being a hasset of peptides with set modifications having that mass
        private readonly Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching;
        private readonly HashSet<DigestionParams> ListOfDigestionParams;

        #endregion Private Fields

        #region Public Constructors

        public ProteinParsimonyEngine(Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, HashSet<DigestionParams> listOfDigestionParams, bool modPeptidesAreDifferent, List<string> nestedIds) : base(nestedIds)
        {
            this.treatModPeptidesAsDifferentPeptides = modPeptidesAreDifferent;
            this.compactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
            this.ListOfDigestionParams = listOfDigestionParams;
        }


        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            ProteinParsimonyResults myAnalysisResults = new ProteinParsimonyResults(this);
            Status("Running protein analysis engine!");

            Status("Applying protein parsimony...");
            myAnalysisResults.ProteinGroups = ApplyProteinParsimony();
            return myAnalysisResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private List<ProteinGroup> ApplyProteinParsimony()
        {
            if (!compactPeptideToProteinPeptideMatching.Values.Any())//if dictionary is empty return an empty list of protein groups
                return new List<ProteinGroup>();

            // digesting an XML database results in a non-mod-agnostic digestion; need to fix this if mod-agnostic parsimony enabled
            if (!treatModPeptidesAsDifferentPeptides)//user want modified and unmodified peptides treated the same
            {
                Dictionary<string, HashSet<PeptideWithSetModifications>> baseSeqToProteinMatch = new Dictionary<string, HashSet<PeptideWithSetModifications>>();
                // dictionary where string key is the base sequence and the HashSet is all PeptidesWithSetModificatiosn with the same sequence 
                // can access which protein these matching peptides came from through the PeptideWithSetModifications object
                foreach (var peptide in compactPeptideToProteinPeptideMatching.SelectMany(b => b.Value))
                {
                    if (baseSeqToProteinMatch.TryGetValue(peptide.BaseSequence, out HashSet<PeptideWithSetModifications> value))
                    {
                        value.Add(peptide);
                    }
                    else
                    {
                        baseSeqToProteinMatch[peptide.BaseSequence] = new HashSet<PeptideWithSetModifications> { peptide };
                    }
                }

                var blah = new Dictionary<PeptideWithSetModifications, List<CompactPeptideBase>>();
                // where to store results
                foreach (var pep in compactPeptideToProteinPeptideMatching)
                {
                    foreach (var pepWithSetMods in pep.Value)
                    {
                        if (blah.TryGetValue(pepWithSetMods, out List<CompactPeptideBase> list))
                        {
                            list.Add(pep.Key);
                        }
                        else
                        {
                            blah.Add(pepWithSetMods, new List<CompactPeptideBase> { pep.Key });
                        }
                    }
                }

                foreach (var baseSequence in baseSeqToProteinMatch)
                {
                    if (baseSequence.Value.Count > 1 && baseSequence.Value.Any(p => p.NumMods > 0))
                    {
                        // list of proteins along with start/end residue in protein and the # missed cleavages
                        var peptideInProteinInfo = new List<Tuple<Protein, DigestionParams, int, int, int>>();
                        foreach (var peptide in baseSequence.Value)// does this peptide have digestion param info?
                            peptideInProteinInfo.Add(new Tuple<Protein, DigestionParams, int, int, int>(peptide.Protein, peptide.digestionParams, peptide.OneBasedStartResidueInProtein, peptide.OneBasedEndResidueInProtein, (int)peptide.MissedCleavages));

                        foreach (var peptide in baseSequence.Value)
                        {
                            foreach (var proteinInfo in peptideInProteinInfo)
                            {

                                var pep = new PeptideWithSetModifications(proteinInfo.Item1, proteinInfo.Item2, proteinInfo.Item3, proteinInfo.Item4, peptide.PeptideDescription, proteinInfo.Item5, peptide.allModsOneIsNterminus, peptide.numFixedMods);
                                foreach (var compactPeptide in blah[peptide])
                                {
                                    compactPeptideToProteinPeptideMatching[compactPeptide].Add(pep);
                                }
                            }
                        }
                    }
                }
            }

            var proteinToPeptidesMatching = new Dictionary<Protein, HashSet<CompactPeptideBase>>();
            var parsimonyProteinList = new Dictionary<Protein, HashSet<CompactPeptideBase>>();
            var proteinsWithUniquePeptides = new Dictionary<Protein, HashSet<PeptideWithSetModifications>>();

            // peptide matched to fullseq (used depending on user preference)
            var compactPeptideToFullSeqMatch = compactPeptideToProteinPeptideMatching.ToDictionary(x => x.Key, x => x.Value.First().Sequence);



            foreach (var kvp in compactPeptideToProteinPeptideMatching)
            {

                HashSet<Protein> proteinsAssociatedWithThisPeptide = new HashSet<Protein>(kvp.Value.Select(p => p.Protein));
                if (proteinsAssociatedWithThisPeptide.Count == 1)
                {
                    if (!proteinsWithUniquePeptides.TryGetValue(kvp.Value.First().Protein, out HashSet<PeptideWithSetModifications> peptides))
                        proteinsWithUniquePeptides.Add(kvp.Value.First().Protein, new HashSet<PeptideWithSetModifications>(kvp.Value));
                    else
                        peptides.UnionWith(kvp.Value);
                }
                //for MultiProtease parsimony protein with the same base sequence can be unique if they 
                //come from different proteolytic digestions
                //Here if proteinAssociatedWithThisPeptide.Count is greater than 1 we must check to ensure that this is not because of
                //different proteolytic digestions yielding the same peptide. 
                //if different proteases do yield the same peptide then both of these peptides are unique
                else
                {

                    foreach (var peptide in kvp.Value)
                    {
                        Protease protease = peptide.digestionParams.Protease;
                        int sameProteaseCount = kvp.Value.Count(v => v.digestionParams.Protease == protease); ;


                        if (sameProteaseCount == 1)
                        {
                            if (!proteinsWithUniquePeptides.TryGetValue(peptide.Protein, out HashSet<PeptideWithSetModifications> peps))
                                proteinsWithUniquePeptides.Add(peptide.Protein, new HashSet<PeptideWithSetModifications> { peptide });
                            else
                                peps.UnionWith(kvp.Value);

                        }
                    }
                }

                // if a peptide is associated with a decoy protein, remove all target protein associations with the peptide
                if (kvp.Value.Any(p => p.Protein.IsDecoy))
                    kvp.Value.RemoveWhere(p => !p.Protein.IsDecoy);

                // if a peptide is associated with a contaminant protein, remove all target protein associations with the peptide
                if (kvp.Value.Any(p => p.Protein.IsContaminant))
                    kvp.Value.RemoveWhere(p => !p.Protein.IsContaminant);

            }


            // makes dictionary with proteins as keys and list of associated peptides as the value (makes parsimony algo easier)
            foreach (var kvp in compactPeptideToProteinPeptideMatching)
            {
                foreach (var peptide in kvp.Value)
                {
                    if (!proteinToPeptidesMatching.TryGetValue(peptide.Protein, out HashSet<CompactPeptideBase> peptides))
                    {
                        proteinToPeptidesMatching.Add(peptide.Protein, new HashSet<CompactPeptideBase>() { kvp.Key });
                    }
                    else
                    {
                        peptides.Add(kvp.Key);
                    }
                }
            }

            // build protein list for each peptide before parsimony has been applied
            var peptideSeqProteinListMatch = new Dictionary<string, HashSet<Protein>>();
            foreach (var kvp in proteinToPeptidesMatching)
            {
                foreach (var peptide in kvp.Value)
                {
                    string pepSequence;
                    if (!treatModPeptidesAsDifferentPeptides)
                    {
                        string nTerminalMasses = peptide.NTerminalMasses == null ? "" : string.Join("", peptide.NTerminalMasses.Select(b => b.ToString(CultureInfo.InvariantCulture)));
                        string cTerminalMasses = peptide.CTerminalMasses == null ? "" : string.Join("", peptide.CTerminalMasses.Select(b => b.ToString(CultureInfo.InvariantCulture)));
                        pepSequence = nTerminalMasses + cTerminalMasses + peptide.MonoisotopicMassIncludingFixedMods.ToString(CultureInfo.InvariantCulture);
                    }
                    else
                    {
                        pepSequence = compactPeptideToFullSeqMatch[peptide];
                    }
                    if (!peptideSeqProteinListMatch.TryGetValue(pepSequence, out HashSet<Protein> proteinListHere))
                    {
                        peptideSeqProteinListMatch.Add(pepSequence, new HashSet<Protein>() { kvp.Key });
                    }
                    else
                    {
                        proteinListHere.Add(kvp.Key);
                    }
                }
            }

            // dictionary associates proteins w/ unused base seqs (list will shrink over time)
            var algDictionary = new Dictionary<Protein, HashSet<string>>();
            foreach (var kvp in peptideSeqProteinListMatch)
            {
                foreach (var protein in kvp.Value)
                {
                    if (algDictionary.TryGetValue(protein, out HashSet<string> newPeptideBaseSeqs))
                    {
                        newPeptideBaseSeqs.Add(kvp.Key);
                    }
                    else
                    {
                        algDictionary.Add(protein, new HashSet<string> { kvp.Key });
                    }
                }
            }

            // dictionary associates proteins w/ unused base seqs (list will NOT shrink over time)
            var proteinToPepSeqMatch = algDictionary.ToDictionary(x => x.Key, x => x.Value);

            // *** main parsimony loop
            bool uniquePeptidesLeft = proteinsWithUniquePeptides.Any();

            int numNewSeqs = algDictionary.Max(p => p.Value.Count);

            while (numNewSeqs != 0)
            {
                var possibleBestProteinList = new List<KeyValuePair<Protein, HashSet<string>>>();

                if (uniquePeptidesLeft)
                {
                    var proteinsWithUniquePeptidesLeft = algDictionary.Where(p => proteinsWithUniquePeptides.ContainsKey(p.Key));
                    if (proteinsWithUniquePeptidesLeft.Any())
                    {
                        possibleBestProteinList.Add(proteinsWithUniquePeptidesLeft.First());
                    }
                    else
                    {
                        uniquePeptidesLeft = false;
                    }
                }

                // gets list of proteins with the most unaccounted-for peptide base sequences
                if (!uniquePeptidesLeft)
                {
                    possibleBestProteinList = algDictionary.Where(p => p.Value.Count == numNewSeqs).ToList();
                }

                Protein bestProtein = possibleBestProteinList.First().Key;
                HashSet<string> newSeqs = new HashSet<string>(algDictionary[bestProtein]);

                // may need to select different protein
                if (possibleBestProteinList.Count > 1)
                {
                    var proteinsWithTheseBaseSeqs = new HashSet<Protein>();

                    foreach (var kvp in possibleBestProteinList)
                    {
                        if (newSeqs.IsSubsetOf(kvp.Value))
                        {
                            proteinsWithTheseBaseSeqs.Add(kvp.Key);
                        }
                    }

                    if (proteinsWithTheseBaseSeqs.Count > 1)
                    {
                        var proteinsOrderedByTotalPeptideCount = new Dictionary<Protein, HashSet<string>>();
                        foreach (var protein in proteinsWithTheseBaseSeqs)
                        {
                            proteinsOrderedByTotalPeptideCount.Add(protein, proteinToPepSeqMatch[protein]);
                        }
                        bestProtein = proteinsOrderedByTotalPeptideCount.OrderByDescending(kvp => kvp.Value.Count).First().Key;
                    }
                }

                parsimonyProteinList.Add(bestProtein, proteinToPeptidesMatching[bestProtein]);

                // remove used peptides from their proteins
                foreach (var newBaseSeq in newSeqs)
                {
                    HashSet<Protein> proteinsWithThisPeptide = peptideSeqProteinListMatch[newBaseSeq];

                    foreach (var protein in proteinsWithThisPeptide)
                    {
                        algDictionary[protein].Remove(newBaseSeq);
                    }
                }
                algDictionary.Remove(bestProtein);
                numNewSeqs = algDictionary.Any() ? algDictionary.Max(p => p.Value.Count) : 0;
            }

            // *** done with parsimony

            // add indistinguishable proteins
            var proteinsGroupedByNumPeptides = proteinToPeptidesMatching.GroupBy(p => p.Value.Count);
            var parsimonyProteinsGroupedByNumPeptides = parsimonyProteinList.GroupBy(p => p.Value.Count);
            var indistinguishableProteins = new ConcurrentDictionary<Protein, HashSet<CompactPeptideBase>>();

            foreach (var group in proteinsGroupedByNumPeptides)
            {
                var parsimonyProteinsWithSameNumPeptides = parsimonyProteinsGroupedByNumPeptides.FirstOrDefault(p => p.Key == group.Key);
                var list = group.ToList();
                if (parsimonyProteinsWithSameNumPeptides != null)
                {
                    Parallel.ForEach(Partitioner.Create(0, list.Count),
                        new ParallelOptions { MaxDegreeOfParallelism = -1 },
                        (range, loopState) =>
                        {
                            for (int i = range.Item1; i < range.Item2; i++)
                            {
                                foreach (var parsimonyProteinWithThisNumPeptides in parsimonyProteinsWithSameNumPeptides)
                                {
                                    if (parsimonyProteinWithThisNumPeptides.Key != list[i].Key
                                    && proteinToPeptidesMatching[parsimonyProteinWithThisNumPeptides.Key].SetEquals(proteinToPeptidesMatching[list[i].Key]))
                                    {
                                        indistinguishableProteins.GetOrAdd(list[i].Key, proteinToPeptidesMatching[list[i].Key]);
                                    }
                                }
                            }
                        }
                    );
                }
            }
            foreach (var protein in indistinguishableProteins)
            {
                parsimonyProteinList.Add(protein.Key, protein.Value);
            }
            //For MultiProtease Parsimony, allows for proteins containing unique peptides for a given protease but not overall
            //to be added back into the protein list
            //For MultiProtease Parsimony!
            //In MultiProtease Parsimony peptides with the same base sequence can be unique if they come from different proteolytic digestions
            //In the formation of parsimonyProteinList duplicate peptides are remove by base sequence, this removes peptides 
            //that are actually unique to a different protein
            //Here we are adding back in these proteins that contain unique peptides to the parsimonyProteinList
            HashSet<Protease> listOfProteases = new HashSet<Protease>();
            foreach (var dp in ListOfDigestionParams)
            {
                if (!listOfProteases.Contains(dp.Protease))
                {
                    listOfProteases.Add(dp.Protease);
                }
            }
            if (listOfProteases.Count > 1)
            {
                HashSet<Protein> parsimonyProteinSet = new HashSet<Protein>();
                foreach (var protein in parsimonyProteinList)
                {
                    parsimonyProteinSet.Add(protein.Key);
                }
                // add back in proteins that contain unique peptides
                foreach (var prot in proteinsWithUniquePeptides)
                {
                    if (!parsimonyProteinSet.Contains(prot.Key))
                    {
                        parsimonyProteinList.Add(prot.Key, proteinToPeptidesMatching[prot.Key]);
                    }

                }

            }
            foreach (var kvp in compactPeptideToProteinPeptideMatching)
            {
                kvp.Value.RemoveWhere(p => !parsimonyProteinList.ContainsKey(p.Protein));
            }

            Status("Finished Parsimony");
            var testing = new HashSet<PeptideWithSetModifications>(compactPeptideToProteinPeptideMatching.Values.SelectMany(p => p));
            return ConstructProteinGroups(new HashSet<PeptideWithSetModifications>(proteinsWithUniquePeptides.Values.SelectMany(p => p)), new HashSet<PeptideWithSetModifications>(compactPeptideToProteinPeptideMatching.Values.SelectMany(p => p)));
        }

        private List<ProteinGroup> ConstructProteinGroups(HashSet<PeptideWithSetModifications> uniquePeptides, HashSet<PeptideWithSetModifications> allPeptides)
        {
            List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
            var proteinToPeptidesMatching = new Dictionary<Protein, HashSet<PeptideWithSetModifications>>();

            foreach (var peptide in allPeptides)
            {
                if (proteinToPeptidesMatching.TryGetValue(peptide.Protein, out HashSet<PeptideWithSetModifications> peptidesHere))
                {
                    peptidesHere.Add(peptide);
                }
                else
                {
                    proteinToPeptidesMatching.Add(peptide.Protein, new HashSet<PeptideWithSetModifications> { peptide });
                }
            }

            // build protein group after parsimony (each group only has 1 protein at this point, indistinguishables will be merged during scoring)
            foreach (var kvp in proteinToPeptidesMatching)
            {
                var allPeptidesHere = proteinToPeptidesMatching[kvp.Key];
                var uniquePeptidesHere = new HashSet<PeptideWithSetModifications>(allPeptidesHere.Where(p => uniquePeptides.Contains(p)));

                proteinGroups.Add(new ProteinGroup(new HashSet<Protein>() { kvp.Key }, allPeptidesHere, uniquePeptidesHere));
            }

            foreach (var proteinGroup in proteinGroups)
            {
                proteinGroup.AllPeptides.RemoveWhere(p => !proteinGroup.Proteins.Contains(p.Protein));
                proteinGroup.DisplayModsOnPeptides = treatModPeptidesAsDifferentPeptides;
            }

            return proteinGroups;
        }

        #endregion Private Methods
    }
}