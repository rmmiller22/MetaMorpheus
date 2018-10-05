using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MassSpectrometry;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.FdrAnalysis
{
    public class FdrAnalysisEngine : MetaMorpheusEngine
    {
        private List<PeptideSpectralMatch> Psms;
        private readonly int MassDiffAcceptorNumNotches;
        private readonly bool UseDeltaScore;
        private readonly double ScoreCutoff;

        public FdrAnalysisEngine(List<PeptideSpectralMatch> psms, int massDiffAcceptorNumNotches, CommonParameters commonParameters, List<string> nestedIds)
            : base(commonParameters, nestedIds)
        {
            Psms = psms;
            MassDiffAcceptorNumNotches = massDiffAcceptorNumNotches;
            UseDeltaScore = commonParameters.UseDeltaScore;
            ScoreCutoff = commonParameters.ScoreCutoff;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            FdrAnalysisResults myAnalysisResults = new FdrAnalysisResults(this);

            Status("Running FDR analysis...");
            DoFalseDiscoveryRateAnalysis(myAnalysisResults);

            myAnalysisResults.PsmsWithin1PercentFdr = Psms.Count(b => b.FdrInfo.QValue < 0.01);

            return myAnalysisResults;
        }

        public static void CalculateEValue(PeptideSpectralMatch psm, MsDataScan scan, CommonParameters commonParameters, Random r,
            int numShuffledDecoys, double precursorMass)
        {
            PeptideWithSetModifications peptide = psm.BestMatchingPeptides.First().Peptide;
            List<double> randomScores = new List<double>();

            foreach (var decoyShuffledPeptide in GetShuffledPeptides(peptide, numShuffledDecoys, r))
            {
                var theorFragments = decoyShuffledPeptide.Fragment(commonParameters.DissociationType, FragmentationTerminus.Both).ToList();
                var matchedFragments = MatchFragmentIons(scan.MassSpectrum, theorFragments, commonParameters, precursorMass);

                double decoyScore = CalculatePeptideScore(scan, matchedFragments, 0);

                randomScores.Add(decoyScore);
            }

            // maximum likelihood estimate is the average of the scores without the target hit
            double maximumLikelihoodEstimate = randomScores.Average();

            double pValue = 1 - SpecialFunctions.GammaLowerRegularized(maximumLikelihoodEstimate, psm.Score);
            
            //double eValue = psm.Counts * (1 - Math.Pow(pValue, psm.Counts));
            double adjustedPValue = 1 - Math.Pow(pValue, psm.Counts);

            // set the PSM's e-value
            psm.SetFdrValues(0, 0, 0, 0, 0, 0, adjustedPValue, -10 * Math.Log(adjustedPValue, 10));
        }

        private void DoFalseDiscoveryRateAnalysis(FdrAnalysisResults myAnalysisResults)
        {
            // Stop if canceled
            if (GlobalVariables.StopLoops) { return; }

            double cumulativeTarget = 0;
            double cumulativeDecoy = 0;

            //Calculate delta scores for the psms (regardless of if we are using them)
            foreach (PeptideSpectralMatch psm in Psms.Where(p => p != null))
            {
                psm.CalculateDeltaScore(ScoreCutoff);
            }

            //determine if Score or DeltaScore performs better
            if (UseDeltaScore)
            {
                const double qValueCutoff = 0.01; //optimize to get the most PSMs at a 1% FDR

                List<PeptideSpectralMatch> scoreSorted = Psms.OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();
                int ScorePSMs = GetNumPSMsAtqValueCutoff(scoreSorted, qValueCutoff);
                scoreSorted = Psms.OrderByDescending(b => b.DeltaScore).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();
                int DeltaScorePSMs = GetNumPSMsAtqValueCutoff(scoreSorted, qValueCutoff);

                //sort by best method
                myAnalysisResults.DeltaScoreImprovement = DeltaScorePSMs > ScorePSMs;
                Psms = myAnalysisResults.DeltaScoreImprovement ?
                    Psms.OrderByDescending(b => b.DeltaScore).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).ToList() :
                    Psms.OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).ToList();
            }
            else //sort by score
            {
                Psms = Psms.OrderByDescending(b => b.Score)
                    .ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).ToList();
            }

            // sort by e-score
            if (commonParameters.CalculateEValue)
            {
                Psms = Psms.OrderByDescending(b => b.FdrInfo.EScore)
                    .ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).ToList();
            }

            //set up arrays for local FDRs
            double[] cumulativeTargetPerNotch = new double[MassDiffAcceptorNumNotches + 1];
            double[] cumulativeDecoyPerNotch = new double[MassDiffAcceptorNumNotches + 1];

            //Assign FDR values to PSMs
            for (int i = 0; i < Psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                var psm = Psms[i];
                int notch = psm.Notch ?? MassDiffAcceptorNumNotches;
                
                //if (psm.IsDecoy)
                //{
                //    cumulativeDecoy++;
                //    cumulativeDecoyPerNotch[notch]++;
                //}
                //else
                //{
                //    cumulativeTarget++;
                //    cumulativeTargetPerNotch[notch]++;
                //}
                cumulativeDecoy += psm.FdrInfo.EValue;

                double qValue = cumulativeDecoy / cumulativeTarget;
                double qValueNotch = cumulativeDecoyPerNotch[notch] / cumulativeTargetPerNotch[notch];

                if (qValue > 1)
                {
                    qValue = 1;
                }
                if (qValueNotch > 1)
                {
                    qValueNotch = 1;
                }

                double eValue = double.NaN;
                double eScore = double.NaN;

                if (psm.FdrInfo != null)
                {
                    eValue = psm.FdrInfo.EValue;
                    eScore = psm.FdrInfo.EScore;
                }

                psm.SetFdrValues((int)cumulativeTarget, (int)cumulativeDecoy, qValue, (int)cumulativeTargetPerNotch[notch], (int)cumulativeDecoyPerNotch[notch], qValueNotch, eValue, eScore);
            }

            //Populate min qValues
            double min_q_value = double.PositiveInfinity;
            double[] min_q_value_notch = new double[MassDiffAcceptorNumNotches + 1];
            for (int i = 0; i < MassDiffAcceptorNumNotches + 1; i++)
            {
                min_q_value_notch[i] = double.PositiveInfinity;
            }

            //The idea here is to set previous qValues as thresholds,
            //such that a lower scoring PSM can't have a higher confidence than a higher scoring PSM
            for (int i = Psms.Count - 1; i >= 0; i--)
            {
                PeptideSpectralMatch psm = Psms[i];
                if (psm.FdrInfo.QValue > min_q_value)
                {
                    psm.FdrInfo.QValue = min_q_value;
                }
                else if (psm.FdrInfo.QValue < min_q_value)
                {
                    min_q_value = psm.FdrInfo.QValue;
                }
                int notch = psm.Notch ?? MassDiffAcceptorNumNotches;
                if (psm.FdrInfo.QValueNotch > min_q_value_notch[notch])
                {
                    psm.FdrInfo.QValueNotch = min_q_value_notch[notch];
                }
                else if (psm.FdrInfo.QValueNotch < min_q_value_notch[notch])
                {
                    min_q_value_notch[notch] = psm.FdrInfo.QValueNotch;
                }
            }
        }

        private static int GetNumPSMsAtqValueCutoff(List<PeptideSpectralMatch> psms, double qValueCutoff)
        {
            int cumulative_target = 0;
            int cumulative_decoy = 0;
            foreach (PeptideSpectralMatch psm in psms)
            {
                if (psm.IsDecoy)
                {
                    cumulative_decoy++;
                    if ((double)cumulative_decoy / cumulative_target >= qValueCutoff)
                    {
                        return cumulative_target;
                    }
                }
                else
                    cumulative_target++;
            }
            return cumulative_target;
        }

        private static IEnumerable<PeptideWithSetModifications> GetShuffledPeptides(PeptideWithSetModifications originalPeptide, int numShuffledPeptidesToGenerate, Random random)
        {
            string[] residues;

            if (originalPeptide.AllModsOneIsNterminus.Any())
            {
                residues = SplitPeptideIntoModifiedResidues(originalPeptide);
            }
            else
            {
                residues = originalPeptide.BaseSequence.Select(p => p.ToString()).ToArray();
            }

            for (int i = 0; i < numShuffledPeptidesToGenerate; i++)
            {
                string shuffledSequence = originalPeptide.FullSequence;

                while (shuffledSequence == originalPeptide.FullSequence && originalPeptide.Length > 3)
                {
                    // Knuth shuffle algorithm
                    // don't shuffle N or C terminal residues
                    for (int t = 1; t < residues.Length - 1; t++)
                    {
                        string tmp = residues[t];
                        int r = random.Next(t, residues.Length - 1);
                        residues[t] = residues[r];
                        residues[r] = tmp;
                    }

                    shuffledSequence = string.Join("", residues);
                }

                yield return new PeptideWithSetModifications(shuffledSequence, GlobalVariables.AllModsKnownDictionary);
            }
        }

        private static string[] SplitPeptideIntoModifiedResidues(PeptideWithSetModifications originalPeptide)
        {
            string[] residues = new string[originalPeptide.BaseSequence.Length];

            for (int i = 0; i < residues.Length; i++)
            {
                residues[i] = "";
            }

            StringBuilder mod = new StringBuilder();
            int aa = 0;
            bool currentlyReadingMod = false;
            int bracketCount = 0;

            for (int r = 0; r < originalPeptide.FullSequence.Length; r++)
            {
                char c = originalPeptide.FullSequence[r];

                switch (c)
                {
                    case '[':
                        currentlyReadingMod = true;
                        mod.Append(c);
                        bracketCount++;
                        break;

                    case ']':
                        bracketCount--;
                        mod.Append(c);

                        if (bracketCount == 0)
                        {
                            currentlyReadingMod = false;

                            if (aa == 0)
                            {
                                residues[aa] += mod.ToString();
                            }
                            else
                            {
                                residues[aa - 1] += mod.ToString();
                            }
                            mod.Clear();
                        }

                        break;

                    default:
                        if (!currentlyReadingMod)
                        {
                            residues[aa] += c.ToString();
                            mod.Clear();
                            aa++;
                        }
                        else
                        {
                            mod.Append(c);
                        }
                        break;
                }
            }

            return residues;
        }
    }
}
