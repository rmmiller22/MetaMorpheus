﻿using MassSpectrometry;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace EngineLayer.Indexing
{
    public class PrecursorIndexingEngine : IndexingEngine
    {
        public PrecursorIndexingEngine(List<Protein> proteinList, List<Modification> variableModifications, List<Modification> fixedModifications, int currentPartition, DecoyType decoyType, IEnumerable<DigestionParams> CollectionOfDigestionParams, CommonParameters commonParams, double maxFragmentSize, List<string> nestedIds) : base(proteinList, variableModifications, fixedModifications, currentPartition, decoyType, CollectionOfDigestionParams, commonParams, maxFragmentSize, nestedIds)
        {
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append("Precursor Mass Only");
            sb.AppendLine("Index partitions: " + CurrentPartition + "/" + commonParameters.TotalPartitions);
            sb.AppendLine("Search Decoys: " + DecoyType);
            sb.AppendLine("Number of proteins: " + ProteinList.Count);
            sb.AppendLine("Number of fixed mods: " + FixedModifications.Count);
            sb.AppendLine("Number of variable mods: " + VariableModifications.Count);
            sb.AppendLine("lp: " + string.Join(",", ProductTypes));
            foreach (var digestionParams in CollectionOfDigestionParams)
            {
                sb.AppendLine("protease: " + digestionParams.Protease);
                sb.AppendLine("initiatorMethionineBehavior: " + digestionParams.InitiatorMethionineBehavior);
                sb.AppendLine("maximumMissedCleavages: " + digestionParams.MaxMissedCleavages);
                sb.AppendLine("minPeptideLength: " + digestionParams.MinPeptideLength);
                sb.AppendLine("maxPeptideLength: " + digestionParams.MaxPeptideLength);
                sb.AppendLine("maximumVariableModificationIsoforms: " + digestionParams.MaxModificationIsoforms);
            }
            sb.Append("Localizeable mods: " + ProteinList.Select(b => b.OneBasedPossibleLocalizedModifications.Count).Sum());
            return sb.ToString();
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            
            // digest database
            HashSet<PeptideWithSetModifications> peptideToId = new HashSet<PeptideWithSetModifications>();

            Parallel.ForEach(Partitioner.Create(0, ProteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, (range, loopState) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops)
                    {
                        loopState.Stop();
                        return;
                    }

                    foreach (var digestionParams in CollectionOfDigestionParams)
                    {
                        foreach (var pepWithSetMods in ProteinList[i].Digest(digestionParams, FixedModifications, VariableModifications))
                        {
                            var observed = peptideToId.Contains(pepWithSetMods);
                            if (observed)
                                continue;
                            lock (peptideToId)
                            {
                                observed = peptideToId.Contains(pepWithSetMods);
                                if (observed)
                                    continue;
                                peptideToId.Add(pepWithSetMods);
                            }
                        }
                    }

                    progress++;
                    var percentProgress = (int)((progress / ProteinList.Count) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Digesting proteins for precursor...", nestedIds));
                    }
                }
            });

            // sort peptides by mass
            var peptidesSortedByMass = peptideToId.AsParallel().WithDegreeOfParallelism(commonParameters.MaxThreadsToUsePerFile).OrderBy(p => p.MonoisotopicMass).ToList();
            peptideToId = null;

            // create fragment index
            int maxFragmentMass = 0;
            for (int i = peptidesSortedByMass.Count - 1; i >= 0; i--)
            {
                if (!Double.IsNaN(peptidesSortedByMass[i].MonoisotopicMass))
                {
                    maxFragmentMass = (int)Math.Min(MaxFragmentSize, (int)Math.Ceiling(Chemistry.ClassExtensions.ToMz(peptidesSortedByMass[i].MonoisotopicMass, 1)));
                    break;
                }
            }

            var fragmentIndex = new List<int>[maxFragmentMass * FragmentBinsPerDalton + 1];

            // populate fragment index
            progress = 0;
            oldPercentProgress = 0;
            for (int i = 0; i < peptidesSortedByMass.Count; i++)
            {
                double mz = Chemistry.ClassExtensions.ToMz(peptidesSortedByMass[i].MonoisotopicMass, 1);
                if (!Double.IsNaN(mz))
                {
                    if (mz < MaxFragmentSize) //if the precursor is larger than the index allows, then stop adding precursors
                    {
                        break;
                    }

                    int fragmentBin = (int)Math.Round(mz * FragmentBinsPerDalton);

                    if (fragmentIndex[fragmentBin] == null)
                        fragmentIndex[fragmentBin] = new List<int> { i };
                    else
                        fragmentIndex[fragmentBin].Add(i);
                }
                progress++;
                var percentProgress = (int)((progress / peptidesSortedByMass.Count) * 100);

                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    ReportProgress(new ProgressEventArgs(percentProgress, "Creating fragment index for precursor...", nestedIds));
                }
            }

            return new IndexingResults(peptidesSortedByMass, fragmentIndex, this);
        }
    }
}