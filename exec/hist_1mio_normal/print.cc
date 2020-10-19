{
   TFile f("hist_no_cut.root");
   TIter next(f.GetListOfKeys());
   TKey *key;
   while ((key=(TKey*)next())) {
      printf("key: %s points to an object of class: %s at %d \n",
      key->GetName(),
      key->GetClassName(),key->GetSeekKey());
   }
}
