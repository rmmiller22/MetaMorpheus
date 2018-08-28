using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public class PepsWSMods: PeptideWithSetModifications
    {
        public readonly Dictionary<int, ModificationWithMass> AllModsOneIsNterminus;
        public readonly DigestionParams DigestionParams;
        public readonly int NumFixedMods;

        //public PepsWSMods(int numFixedMods, Protein protein, DigestionParams digestionParams,int proteinOneBasedStart, int proteinOneBasedEnd, Dictionary<int, ModificationWithMass> allModsOneIsNterminus = null, int missedCleavages = 0)
        //    : base(numFixedMods, protein, proteinOneBasedStart, proteinOneBasedEnd, allModsOneIsNterminus, missedCleavages)
        //{
        //    NumFixedMods = numFixedMods;
        //    AllModsOneIsNterminus = allModsOneIsNterminus ?? new Dictionary<int, ModificationWithMass>();
        //    DigestionParams = digestionParams;
        //}

            //test
        public PepsWSMods(Protein protein, DigestionParams digestionParams, int oneBasedStartResidueInProtein, int oneBasedEndResidueInProtein, string peptideDescription, int missedCleavages, Dictionary<int, ModificationWithMass> allModsOneIsNterminus, int numFixedMods)
            : base(protein, digestionParams, oneBasedStartResidueInProtein, oneBasedEndResidueInProtein, peptideDescription, missedCleavages, allModsOneIsNterminus, numFixedMods)
        {
            AllModsOneIsNterminus = allModsOneIsNterminus;
            NumFixedMods = numFixedMods;
            DigestionParams = digestionParams;

        }

        public override bool Equals(object obj)
        {
            var q = obj as PeptideWithSetModifications;
            return q != null
                && q.Sequence.Equals(Sequence)
                && q.OneBasedEndResidueInProtein == OneBasedEndResidueInProtein
                && (q.Protein.Accession == null && Protein.Accession == null) || q.Protein.Accession.Equals(Protein.Accession)
                && q.DigestionParams.Protease == DigestionParams.Protease;
        }

        public override int GetHashCode()
        {
            return Sequence.GetHashCode() + DigestionParams.Protease.GetHashCode();
  
        }


    }
}
