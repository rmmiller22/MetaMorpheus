using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class BinGenerationTest
    {
        [Test]
        public static void TestBinGeneration()
        {
            SearchTask st = new SearchTask
            {
                CommonParameters = new CommonParameters(scoreCutoff: 1, digestionParams: new DigestionParams(minPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain)),

                SearchParameters = new SearchParameters
                {
                    DoHistogramAnalysis = true,
                    MassDiffAcceptorType = MassDiffAcceptorType.Open,
                    DecoyType = DecoyType.None,
                    DoParsimony = true,
                    DoQuantification = true
                },
            };

            string proteinDbFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "BinGenerationTest.xml");
            string mzmlFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "BinGenerationTest.mzML");

            Protein prot1 = new Protein("MEDEEK", "prot1");
            Protein prot2 = new Protein("MENEEK", "prot2");

            ModificationMotif.TryGetMotif("D", out ModificationMotif motif);
            ModificationWithMass mod = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, 10);

            var pep1_0 = prot1.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).First();
            PepsWSMods pwsm1_0 = new PepsWSMods(pep1_0.Protein, pep1_0.DigestionParams, pep1_0.OneBasedStartResidueInProtein, pep1_0.OneBasedEndResidueInProtein, pep1_0.PeptideDescription, pep1_0.MissedCleavages, pep1_0.AllModsOneIsNterminus, pep1_0.NumFixedMods);
            var pep1_10 = prot1.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).Last();
            PepsWSMods pwsm1_10 = new PepsWSMods(pep1_10.Protein, pep1_10.DigestionParams, pep1_10.OneBasedStartResidueInProtein, pep1_10.OneBasedEndResidueInProtein, pep1_10.PeptideDescription, pep1_10.MissedCleavages, pep1_10.AllModsOneIsNterminus, pep1_10.NumFixedMods);

            Protein prot3 = new Protein("MAAADAAAAAAAAAAAAAAA", "prot3");

            var pep2_0 = prot3.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).First();
            PepsWSMods pwsm2_0 = new PepsWSMods(pep2_0.Protein, pep2_0.DigestionParams, pep2_0.OneBasedStartResidueInProtein, pep2_0.OneBasedEndResidueInProtein, pep2_0.PeptideDescription, pep2_0.MissedCleavages, pep2_0.AllModsOneIsNterminus, pep2_0.NumFixedMods);
            var pep2_10 = prot3.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass> { mod }).Last();
            PepsWSMods pwsm2_10 = new PepsWSMods(pep2_10.Protein, pep2_10.DigestionParams, pep2_10.OneBasedStartResidueInProtein, pep2_10.OneBasedEndResidueInProtein, pep2_10.PeptideDescription, pep2_10.MissedCleavages, pep2_10.AllModsOneIsNterminus, pep2_10.NumFixedMods);

            Protein prot4 = new Protein("MNNDNNNN", "prot4");
            var pep3_10 = prot4.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass> { mod }).Last();
            PepsWSMods pwsm3_10 = new PepsWSMods(pep3_10.Protein, pep3_10.DigestionParams, pep3_10.OneBasedStartResidueInProtein, pep3_10.OneBasedEndResidueInProtein, pep3_10.PeptideDescription, pep3_10.MissedCleavages, pep3_10.AllModsOneIsNterminus, pep3_10.NumFixedMods);

            List<PepsWSMods> pepsWithSetMods = new List<PepsWSMods> { pwsm1_0, pwsm1_10, pwsm2_0, pwsm2_10, pwsm3_10 };
            MsDataFile myMsDataFile = new TestDataFile(pepsWithSetMods);

            List<Protein> proteinList = new List<Protein> { prot1, prot2, prot3, prot4 };

            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlFilePath, false);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinList, proteinDbFilePath);

            string output_folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestBinGeneration");
            Directory.CreateDirectory(output_folder);
            st.RunTask(
                output_folder,
                new List<DbForTask> { new DbForTask(proteinDbFilePath, false) },
                new List<string> { mzmlFilePath },
                null);

            Assert.AreEqual(3, File.ReadLines(Path.Combine(output_folder, @"MassDifferenceHistogram.tsv")).Count());
        }

        [Test]
        public static void TestProteinSplitAcrossFiles()
        {
            SearchTask st = new SearchTask()
            {
                CommonParameters = new CommonParameters(
                    scoreCutoff: 1,
                    digestionParams: new DigestionParams(
                        maxMissedCleavages: 0,
                        minPeptideLength: 5,
                        initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain)),

                SearchParameters = new SearchParameters
                {
                    DoHistogramAnalysis = true,
                    MassDiffAcceptorType = MassDiffAcceptorType.Open,
                    MatchBetweenRuns = true,
                    DoQuantification = true
                },
            };

            string proteinDbFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestProteinSplitAcrossFiles.xml");
            string mzmlFilePath1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestProteinSplitAcrossFiles1.mzML");
            string mzmlFilePath2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestProteinSplitAcrossFiles2.mzML");

            ModificationMotif.TryGetMotif("D", out ModificationMotif motif);
            ModificationWithMass mod = new ModificationWithMass("mod1", "mt", motif, TerminusLocalization.Any, 10);

            IDictionary<int, List<Modification>> oneBasedModification = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification>{ mod } }
            };

            Protein prot1 = new Protein("MEDEEK", "prot1", oneBasedModifications: oneBasedModification);

            var pep1 = prot1.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).First();
            PepsWSMods pwsm1 = new PepsWSMods(pep1.Protein, pep1.DigestionParams, pep1.OneBasedStartResidueInProtein, pep1.OneBasedEndResidueInProtein, pep1.PeptideDescription, pep1.MissedCleavages, pep1.AllModsOneIsNterminus, pep1.NumFixedMods);
            var pep2 = prot1.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).Last();
            PepsWSMods pwsm2 = new PepsWSMods(pep1.Protein, pep1.DigestionParams, pep1.OneBasedStartResidueInProtein, pep1.OneBasedEndResidueInProtein, pep1.PeptideDescription, pep1.MissedCleavages, pep1.AllModsOneIsNterminus, pep1.NumFixedMods);

            List<PepsWSMods> listForFile1 = new List<PepsWSMods> { pwsm1, pwsm2 };
            List<PepsWSMods> listForFile2 = new List<PepsWSMods> { pwsm2 };
            MsDataFile myMsDataFile1 = new TestDataFile(listForFile1);
            MsDataFile myMsDataFile2 = new TestDataFile(listForFile2);

            List<Protein> proteinList = new List<Protein> { prot1 };

            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlFilePath1, false);
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile2, mzmlFilePath2, false);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinList, proteinDbFilePath);

            string output_folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestProteinSplitAcrossFiles");
            Directory.CreateDirectory(output_folder);

            st.RunTask(
                output_folder,
                new List<DbForTask> { new DbForTask(proteinDbFilePath, false) },
                new List<string> { mzmlFilePath1, mzmlFilePath2, },
                null);
        }
    }
}