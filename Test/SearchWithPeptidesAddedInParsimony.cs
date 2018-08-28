﻿using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class SearchWithPeptidesAddedInParsimony
    {
        [Test]
        public static void SearchWithPeptidesAddedInParsimonyTest()
        {
            // Make sure can run the complete search task when multiple compact peptides may correspond to a single PWSM
            SearchTask st = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    DecoyType = DecoyType.None,
                    ModPeptidesAreDifferent = false
                },
                CommonParameters = new CommonParameters(scoreCutoff: 1, digestionParams: new DigestionParams(minPeptideLength: 2)),
            };

            string xmlName = "andguiaheow.xml";

            CommonParameters CommonParameters = new CommonParameters(
                scoreCutoff: 1,
                digestionParams: new DigestionParams(
                    maxMissedCleavages: 0,
                    minPeptideLength: 1,
                    maxModificationIsoforms: 2,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                    maxModsForPeptides: 1));

            ModificationMotif.TryGetMotif("A", out ModificationMotif motifA);
            ModificationWithMass alanineMod = new ModificationWithMass("111", "mt", motifA, TerminusLocalization.Any, 111);

            var variableModifications = new List<ModificationWithMass>();
            IDictionary<int, List<Modification>> oneBasedModifications1 = new Dictionary<int, List<Modification>>
                {
                    {2, new List<Modification>{ alanineMod } }
                };
            Protein protein1 = new Protein("MA", "protein1", oneBasedModifications: oneBasedModifications1);
            // Alanine = Glycine + CH2

            ModificationMotif.TryGetMotif("G", out ModificationMotif motif1);

            ModificationWithMass glycineMod = new ModificationWithMass("CH2 on Glycine", "mt", motif1, TerminusLocalization.Any, Chemistry.ChemicalFormula.ParseFormula("CH2").MonoisotopicMass);

            IDictionary<int, List<Modification>> oneBasedModifications2 = new Dictionary<int, List<Modification>>
                {
                    {2, new List<Modification>{glycineMod} }
                };
            Protein protein2 = new Protein("MG", "protein3", oneBasedModifications: oneBasedModifications2);

            PeptideWithSetModifications pwsmMA = protein1.Digest(CommonParameters.DigestionParams, new List<ModificationWithMass>(), variableModifications).First();
            PepsWSMods pepMA = new PepsWSMods(pwsmMA.Protein, pwsmMA.DigestionParams, pwsmMA.OneBasedStartResidueInProtein, pwsmMA.OneBasedEndResidueInProtein, pwsmMA.PeptideDescription, pwsmMA.MissedCleavages, pwsmMA.AllModsOneIsNterminus, pwsmMA.NumFixedMods);
            PeptideWithSetModifications pwsmMA111 = protein1.Digest(CommonParameters.DigestionParams, new List<ModificationWithMass>(), variableModifications).Last();
            PepsWSMods pepMA111 = new PepsWSMods(pwsmMA111.Protein, pwsmMA111.DigestionParams, pwsmMA111.OneBasedStartResidueInProtein, pwsmMA111.OneBasedEndResidueInProtein, pwsmMA111.PeptideDescription, pwsmMA111.MissedCleavages, pwsmMA111.AllModsOneIsNterminus, pwsmMA111.NumFixedMods);
            var pwsmMG = protein2.Digest(CommonParameters.DigestionParams, new List<ModificationWithMass>(), variableModifications).First();
            PepsWSMods pepMG = new PepsWSMods(pwsmMG.Protein, pwsmMG.DigestionParams, pwsmMG.OneBasedStartResidueInProtein, pwsmMG.OneBasedEndResidueInProtein, pwsmMG.PeptideDescription, pwsmMG.MissedCleavages, pwsmMG.AllModsOneIsNterminus, pwsmMG.NumFixedMods);

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { protein1, protein2 }, xmlName);

            string mzmlName = @"ajgdiu.mzML";

            MsDataFile myMsDataFile = new TestDataFile(new List<PepsWSMods> { pepMA, pepMG, pepMA111 }, true);

            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

            st.RunTask("",
                new List<DbForTask> { new DbForTask(xmlName, false) },
                new List<string> { mzmlName }, "");
        }
    }
}