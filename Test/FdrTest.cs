using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class FdrTest
    {
        [Test]
        public static void FdrTestMethod()
        {
            MassDiffAcceptor searchModes = new DotMassDiffAcceptor(null, new List<double> { 0, 1.0029 }, new PpmTolerance(5));
            List<string> nestedIds = new List<string>();

            Protein p = new Protein("MNKNNKNNNKNNNNK", null);
            DigestionParams digestionParams = new DigestionParams();
            var digested = p.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();

            PepsWSMods pep1 = new PepsWSMods(digested[0].Protein, digested[0].DigestionParams, digested[0].OneBasedStartResidueInProtein, digested[0].OneBasedEndResidueInProtein, digested[0].PeptideDescription, digested[0].MissedCleavages, digested[0].AllModsOneIsNterminus, digested[0].NumFixedMods);
            PepsWSMods pep2 = new PepsWSMods(digested[1].Protein, digested[1].DigestionParams, digested[1].OneBasedStartResidueInProtein, digested[1].OneBasedEndResidueInProtein, digested[1].PeptideDescription, digested[1].MissedCleavages, digested[1].AllModsOneIsNterminus, digested[1].NumFixedMods);
            PepsWSMods pep3 = new PepsWSMods(digested[2].Protein, digested[2].DigestionParams, digested[2].OneBasedStartResidueInProtein, digested[2].OneBasedEndResidueInProtein, digested[2].PeptideDescription, digested[2].MissedCleavages, digested[2].AllModsOneIsNterminus, digested[2].NumFixedMods);
            PepsWSMods pep4 = new PepsWSMods(digested[3].Protein, digested[3].DigestionParams, digested[3].OneBasedStartResidueInProtein, digested[3].OneBasedEndResidueInProtein, digested[3].PeptideDescription, digested[3].MissedCleavages, digested[3].AllModsOneIsNterminus, digested[3].NumFixedMods);
            
            TestDataFile t = new TestDataFile(new List<PepsWSMods> { pep1, pep2, pep3 });

            CompactPeptide peptide1 = new CompactPeptide(pep1, TerminusType.None);
            MsDataScan mzLibScan1 = t.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scan1 = new Ms2ScanWithSpecificMass(mzLibScan1, peptide1.MonoisotopicMassIncludingFixedMods.ToMz(1), 1, null);
            PeptideSpectralMatch psm1 = new PeptideSpectralMatch(peptide1, 0, 3, 0, scan1, digestionParams);

            CompactPeptide peptide2 = new CompactPeptide(pep2, TerminusType.None);
            MsDataScan mzLibScan2 = t.GetOneBasedScan(4);
            Ms2ScanWithSpecificMass scan2 = new Ms2ScanWithSpecificMass(mzLibScan2, peptide2.MonoisotopicMassIncludingFixedMods.ToMz(1), 1, null);
            PeptideSpectralMatch psm2 = new PeptideSpectralMatch(peptide2, 1, 2, 1, scan2, digestionParams);

            CompactPeptide peptide3 = new CompactPeptide(pep3, TerminusType.None);
            MsDataScan mzLibScan3 = t.GetOneBasedScan(6);
            Ms2ScanWithSpecificMass scan3 = new Ms2ScanWithSpecificMass(mzLibScan3, peptide3.MonoisotopicMassIncludingFixedMods.ToMz(1), 1, null);
            PeptideSpectralMatch psm3 = new PeptideSpectralMatch(peptide3, 0, 1, 2, scan3, digestionParams);

            CompactPeptide peptide4 = new CompactPeptide(pep4, TerminusType.None);
            psm3.AddOrReplace(peptide4, 1, 1, true);

            Dictionary<CompactPeptideBase, HashSet<PepsWSMods>> matching = new Dictionary<CompactPeptideBase, HashSet<PepsWSMods>>
            {
                {
                    peptide1, new HashSet<PepsWSMods>{ pep1 }
                },
                {
                    peptide2, new HashSet<PepsWSMods>{ pep2 }
                },
                {
                    peptide3, new HashSet<PepsWSMods>{ pep3 }
                },
                {
                    peptide4, new HashSet<PepsWSMods>{ pep4 }
                },
            };

            psm1.MatchToProteinLinkedPeptides(matching);
            psm2.MatchToProteinLinkedPeptides(matching);
            psm3.MatchToProteinLinkedPeptides(matching);

            var newPsms = new List<PeptideSpectralMatch> { psm1, psm2, psm3 };
            CommonParameters cp = new CommonParameters(calculateEValue: true);

            FdrAnalysisEngine fdr = new FdrAnalysisEngine(newPsms, searchModes.NumNotches, cp, nestedIds);

            fdr.Run();

            Assert.AreEqual(2, searchModes.NumNotches);
            Assert.AreEqual(0, newPsms[0].FdrInfo.CumulativeDecoyNotch);
            Assert.AreEqual(1, newPsms[0].FdrInfo.CumulativeTargetNotch);
            Assert.AreEqual(0, newPsms[1].FdrInfo.CumulativeDecoyNotch);
            Assert.AreEqual(1, newPsms[1].FdrInfo.CumulativeTargetNotch);
            Assert.AreEqual(0, newPsms[2].FdrInfo.CumulativeDecoyNotch);
            Assert.AreEqual(1, newPsms[2].FdrInfo.CumulativeTargetNotch);

            Assert.AreEqual(0, newPsms[0].FdrInfo.CumulativeDecoy);
            Assert.AreEqual(1, newPsms[0].FdrInfo.CumulativeTarget);
            Assert.AreEqual(0, newPsms[1].FdrInfo.CumulativeDecoy);
            Assert.AreEqual(2, newPsms[1].FdrInfo.CumulativeTarget);
            Assert.AreEqual(0, newPsms[2].FdrInfo.CumulativeDecoy);
            Assert.AreEqual(3, newPsms[2].FdrInfo.CumulativeTarget);
        }

        [Test]
        public static void TestDeltaValues()
        {
            CommonParameters CommonParameters = new CommonParameters(scoreCutoff: 1, useDeltaScore: true, digestionParams: new DigestionParams(minPeptideLength: 5));

            SearchParameters SearchParameters = new SearchParameters
            {
                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
            };
            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();

            // Generate data for files
            Protein TargetProtein1 = new Protein("TIDEANTHE", "accession1");
            Protein TargetProtein2 = new Protein("TIDELVE", "accession2");
            Protein TargetProtein3 = new Protein("TIDENIE", "accession3");
            Protein TargetProteinLost = new Protein("PEPTIDEANTHE", "accession4");
            Protein DecoyProteinFound = new Protein("PETPLEDQGTHE", "accessiond", isDecoy: true);
            var pep1 = TargetProtein1.Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0];
            PepsWSMods pwsm1 = new PepsWSMods(pep1.Protein, pep1.DigestionParams, pep1.OneBasedStartResidueInProtein, pep1.OneBasedEndResidueInProtein, pep1.PeptideDescription, pep1.MissedCleavages, pep1.AllModsOneIsNterminus, pep1.NumFixedMods);

            var pep2 = TargetProtein2.Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0];
            PepsWSMods pwsm2 = new PepsWSMods(pep2.Protein, pep2.DigestionParams, pep2.OneBasedStartResidueInProtein, pep2.OneBasedEndResidueInProtein, pep2.PeptideDescription, pep2.MissedCleavages, pep2.AllModsOneIsNterminus, pep2.NumFixedMods);

            var pep3 = TargetProtein3.Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0];
            PepsWSMods pwsm3 = new PepsWSMods(pep3.Protein, pep3.DigestionParams, pep3.OneBasedStartResidueInProtein, pep3.OneBasedEndResidueInProtein, pep3.PeptideDescription, pep3.MissedCleavages, pep3.AllModsOneIsNterminus, pep3.NumFixedMods);

            var pep4 = DecoyProteinFound.Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0];
            PepsWSMods pwsm4 = new PepsWSMods(pep4.Protein, pep4.DigestionParams, pep4.OneBasedStartResidueInProtein, pep4.OneBasedEndResidueInProtein, pep4.PeptideDescription, pep4.MissedCleavages, pep4.AllModsOneIsNterminus, pep4.NumFixedMods);

            MsDataFile myMsDataFile = new TestDataFile(new List<PepsWSMods>
            {
                pwsm1,
                pwsm2,
                pwsm3,
                pwsm4
            });

            var proteinList = new List<Protein> { TargetProtein1, TargetProtein2, TargetProtein3, TargetProteinLost, DecoyProteinFound };

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            //check better when using delta
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, CommonParameters, new List<string>()).Run();

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType>
            { ProductType.B, ProductType.Y }, 1, DecoyType.None, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, 30000, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();
            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);
            PeptideSpectralMatch[] allPsmsArrayModern = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArrayModern, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, massDiffAcceptor, 0, new List<string>()).Run();

            Dictionary<CompactPeptideBase, HashSet<PepsWSMods>> compactPeptideToProteinPeptideMatching =
                new Dictionary<CompactPeptideBase, HashSet<PepsWSMods>>();
            if (proteinList.Any())
            {
                SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(allPsmsArray.ToList(), proteinList, fixedModifications, variableModifications, new List<ProductType>
                { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>());
                var res = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
                compactPeptideToProteinPeptideMatching = res.CompactPeptideToProteinPeptideMatching;
            }

            foreach (var psm in allPsmsArray)
            {
                psm.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);
            }
            foreach (var psm in allPsmsArrayModern)
            {
                psm.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);
            }
            FdrAnalysisResults fdrResultsClassicDelta = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArray.ToList(), 1, CommonParameters, new List<string>()).Run());
            FdrAnalysisResults fdrResultsModernDelta = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArrayModern.ToList(), 1, CommonParameters, new List<string>()).Run());
            Assert.IsTrue(fdrResultsClassicDelta.PsmsWithin1PercentFdr == 3);
            Assert.IsTrue(fdrResultsModernDelta.PsmsWithin1PercentFdr == 3);

            CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5));

            //check worse when using score
            FdrAnalysisResults fdrResultsClassic = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArray.ToList(), 1, CommonParameters, new List<string>()).Run());
            FdrAnalysisResults fdrResultsModern = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArray.ToList(), 1, CommonParameters, new List<string>()).Run());
            Assert.IsTrue(fdrResultsClassic.PsmsWithin1PercentFdr == 0);
            Assert.IsTrue(fdrResultsModern.PsmsWithin1PercentFdr == 0);

            //check that when delta is bad, we used the score
            // Generate data for files
            Protein DecoyProtein1 = new Protein("TLEDAGGTHE", "accession1d", isDecoy: true);
            Protein DecoyProtein2 = new Protein("TLEDLVE", "accession2d", isDecoy: true);
            Protein DecoyProtein3 = new Protein("TLEDNIE", "accession3d", isDecoy: true);
            Protein DecoyProteinShiny = new Protein("GGGGGG", "accessionShinyd", isDecoy: true);
                      

            var pep5 = DecoyProteinShiny.Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0];
            PepsWSMods pwsm5 = new PepsWSMods(pep5.Protein, pep5.DigestionParams, pep5.OneBasedStartResidueInProtein, pep5.OneBasedEndResidueInProtein, pep5.PeptideDescription, pep5.MissedCleavages, pep5.AllModsOneIsNterminus, pep5.NumFixedMods);


            myMsDataFile = new TestDataFile(new List<PepsWSMods>
            {
                pwsm1,
                pwsm2,
                pwsm3,
                pwsm5
            });

            proteinList = new List<Protein>
            {
                TargetProtein1, DecoyProtein1,
                TargetProtein2, DecoyProtein2,
                TargetProtein3, DecoyProtein3,
                DecoyProteinShiny,
            };

            listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            //check no change when using delta
            allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, CommonParameters, new List<string>()).Run();

            CommonParameters = new CommonParameters(useDeltaScore: true, digestionParams: new DigestionParams(minPeptideLength: 5));

            indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType>
            { ProductType.B, ProductType.Y }, 1, DecoyType.None, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, 30000, new List<string>());
            indexResults = (IndexingResults)indexEngine.Run();
            massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);
            allPsmsArrayModern = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArrayModern, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, massDiffAcceptor, 0, new List<string>()).Run();

            var compactPeptideToProteinPeptideMatching2 = new Dictionary<CompactPeptideBase, HashSet<PepsWSMods>>();
            if (proteinList.Any())
            {
                SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine2 = new SequencesToActualProteinPeptidesEngine(allPsmsArray.ToList(), proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>());
                var res = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine2.Run();
                compactPeptideToProteinPeptideMatching2 = res.CompactPeptideToProteinPeptideMatching;
            }

            foreach (var psm in allPsmsArray)
            {
                psm.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching2);
            }
            foreach (var psm in allPsmsArrayModern)
            {
                psm.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching2);
            }
            fdrResultsClassicDelta = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArray.ToList(), 1, CommonParameters, new List<string>()).Run());
            fdrResultsModernDelta = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArrayModern.ToList(), 1, CommonParameters, new List<string>()).Run());
            Assert.IsTrue(fdrResultsClassicDelta.PsmsWithin1PercentFdr == 3);
            Assert.IsTrue(fdrResultsModernDelta.PsmsWithin1PercentFdr == 3);

            CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5));

            //check no change when using score
            fdrResultsClassic = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArray.ToList(), 1, CommonParameters, new List<string>()).Run());
            fdrResultsModern = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArrayModern.ToList(), 1, CommonParameters, new List<string>()).Run());
            Assert.IsTrue(fdrResultsClassic.PsmsWithin1PercentFdr == 3);
            Assert.IsTrue(fdrResultsModern.PsmsWithin1PercentFdr == 3);
        }
    }
}