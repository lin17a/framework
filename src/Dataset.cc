// -*- C++ -*-
// author: afiq anuar
// short: please refer to header for information

template <typename Tree>
Framework::Dataset<Tree>::Dataset(const std::string &name_, const std::string &tree_name_, 
                                  const std::string &tree_struct_ /*= ""*/, char tree_delim_ /*= ' '*/,
                                  const std::vector<std::string> &v_file_ /*= {}*/) :
  name(name_),
  tree_name(tree_name_),
  tree_struct(tree_struct_),
  tree_delim(tree_delim_),
  v_file(v_file_)
{
  tree_ptr = nullptr;
  v_weight = {};

  if (!v_file.empty())
    evaluate();
}



template <typename Tree>
Framework::Dataset<Tree>::~Dataset()
{
  reset();
}



template <typename Tree>
void Framework::Dataset<Tree>::set_tree(const std::string &tree_name_)
{
  tree_name = tree_name_;
}



template <typename Tree>
void Framework::Dataset<Tree>::set_structure(const std::string &tree_struct_)
{
  tree_struct = tree_struct_;
}



template <typename Tree>
void Framework::Dataset<Tree>::set_delimiter(char tree_delim_)
{
  tree_delim = tree_delim_;
}



template <typename Tree>
void Framework::Dataset<Tree>::set_files(const std::vector<std::string> &v_file_, int nfile /*= -1*/, bool force_replace /*= false*/)
{
  if (v_file_.empty())
    return;

  if (force_replace)
    reset();

  if (!v_file.empty())
    return;

  nfile = (nfile > 0 and nfile <= v_file_.size()) ? nfile : v_file_.size();
  if (nfile != v_file_.size()) {
    for (int ifile = 0; ifile < nfile; ++ifile)
      v_file.emplace_back(v_file_[ifile]);
  }
  else 
    v_file = v_file_;

  evaluate();
}



template <typename Tree>
bool Framework::Dataset<Tree>::add_file(const std::string &file)
{
  auto iF = std::find(std::begin(v_file), std::end(v_file), file);
  if (iF != std::end(v_file))
    return false;

  v_file.emplace_back(file);
  evaluate(-1);
  return true;
}



template <typename Tree>
int Framework::Dataset<Tree>::n_files()
{
  return v_file.size();
}



template <typename Tree>
bool Framework::Dataset<Tree>::add_weight(const std::string &wgt_name, double wgt)
{
  auto iW = std::find_if(std::begin(v_weight), std::end(v_weight), [&wgt_name] (const auto &weight) {return weight.first == wgt_name;});
  if (iW != std::end(v_weight))
    return false;

  v_weight.emplace_back(wgt_name, wgt);
  return true;
}



template <typename Tree>
Framework::Dataset<Tree> Framework::Dataset<Tree>::split()
{
  std::vector<std::string> vf1, vf2;
  const int half = (v_file.size() % 2 == 1) ? int(v_file.size() / 2) + 1 : v_file.size() / 2;

  for (int iF = 0; iF < v_file.size(); ++iF) {
    if (iF <= half)
      vf1.emplace_back(v_file[iF]);
    else
      vf2.emplace_back(v_file[iF]);
  }

  auto dname = name;
  auto vw = v_weight;
  this->name = dname + "_1";
  this->set_files(vf1, true);
  this->v_weight = vw;

  auto dat = Dataset(dname + "_2", tree_name, tree_struct, tree_delim, vf2);
  dat->v_weight = vw;

  return dat;
}



template <typename Tree>
long long Framework::Dataset<Tree>::current_entry(long long entry) const
{
  if (tree_ptr == nullptr)
    return -1;

  return tree_ptr->LoadTree(entry);
}



template <typename Tree>
const std::unique_ptr<Tree>& Framework::Dataset<Tree>::tree() const
{
  return tree_ptr;
}



template <typename Tree>
double Framework::Dataset<Tree>::get_weight(const std::string &wgt_name) const
{
  auto iW = std::find_if(std::begin(v_weight), std::end(v_weight), [&wgt_name] (const auto &weight) {return weight.first == wgt_name;});
  if (iW == std::end(v_weight)) {
    // TODO something about logger telling absence of weight
    return 0.;
  }
  else 
    return iW->second;
}



template <>
void Framework::Dataset<TChain>::evaluate(int index /*= 0*/)
{
  if (tree_ptr == nullptr) {
    tree_ptr = std::make_unique<TChain>(tree_name.c_str());
    //tree_ptr->SetImplicitMT(false);
    tree_ptr->SetBranchStatus("*", 0);
  }

  if (index == 0) {
    for (const auto &file : v_file)
      tree_ptr->Add(file.c_str());
  }
  else if (index == -1)
    tree_ptr->Add(v_file.back().c_str());
}



template <>
void Framework::Dataset<TTree>::evaluate(int index /*= 0*/)
{
  static bool flag_struct = false;
  if (!flag_struct and tree_struct == "") {
    // TODO something about logging the error
    return;
  }

  if (tree_ptr == nullptr) {
    tree_ptr = std::make_unique<TTree>(tree_name.c_str(), "");
    //tree_ptr->SetImplicitMT(false);
    tree_ptr->SetBranchStatus("*", 0);
  }

  if (index == 0) {
    for (const auto &file : v_file) {
      if (!flag_struct) {
        tree_ptr->ReadFile(file.c_str(), tree_struct.c_str(), tree_delim);
        flag_struct = true;
      }
      else
        tree_ptr->ReadFile(file.c_str());
    }
  }
  else if (index == -1) {
    if (!flag_struct) {
      tree_ptr->ReadFile(v_file.back().c_str(), tree_struct.c_str(), tree_delim);
      flag_struct = true;
    }
    else
      tree_ptr->ReadFile(v_file.back().c_str());
  }
}



template <typename Tree>
template <typename ...Collections>
void Framework::Dataset<Tree>::associate(Collections &...colls)
{
  static_assert(sizeof...(colls) > 0, "ERROR: Dataset::allocate makes no sense when called without arguments!!");

  if (tree_ptr == nullptr)
    throw std::runtime_error( "ERROR: Dataset::associate should not be called before assigning the files to be analyzed!!" );

  // apparently the tree is loaded only after this call
  // so associate will fail without it
  tree_ptr->GetEntries();

  (colls.associate(*this), ...);
  allocator.set_allocator([&colls...] () { (colls.reassociate(), ...); });
  tree_ptr->SetNotify(&allocator);
}



template <typename Tree>
template <typename Analyzer>
void Framework::Dataset<Tree>::set_analyzer(Analyzer analyzer_)
{
  using Traits = function_traits<decltype(analyzer_)>;
  static_assert(Traits::arity == 1 and std::is_convertible_v<typename Traits::template bare_arg<0>, long long>, 
                "ERROR: Dataset::set_analyzer only takes one argument that is convertible to entry number!!");
  static_assert(std::is_same_v<typename Traits::result_type, void>, 
                "ERROR: Dataset::set_analyzer: currently non-void return type is not supported!!");

  if (!analyzer)
    analyzer = std::function<void(long long)>(analyzer_);
}



template <typename Tree>
void Framework::Dataset<Tree>::analyze(long long total /*= -1LL*/, long long skip /*= -1LL*/) const
{
  if (tree_ptr == nullptr)
    throw std::runtime_error( "ERROR: Dataset::analyze should not be called before assigning the files to be analyzed!!" );

  if (!allocator)
    throw std::runtime_error( "ERROR: Dataset::analyze should not be called before calling Dataset::associate!!" );

  if (!analyzer)
    throw std::runtime_error( "ERROR: Dataset::analyze should not be called before calling Dataset::set_analyzer!!" );

  const auto dEvt = (total > 0LL and total <= tree_ptr->GetEntries()) ? total : tree_ptr->GetEntries();
  std::cout << "Processing " << dEvt << " events..." << std::endl;

  const bool doskip = skip > 0LL and skip < total;
  if (doskip)
    std::cout << "Skipping " << skip << " events..." << std::endl;

  for (auto cEvt = (doskip) ? skip : 0LL; cEvt < dEvt; ++cEvt)
    analyzer(current_entry(cEvt));
  std::cout << "Processed " << dEvt << " events!" << std::endl;
}



template <typename Tree>
void Framework::Dataset<Tree>::reset()
{
  tree_ptr->ResetBranchAddresses();
  tree_ptr->Reset();

  v_file.clear();
  v_file.shrink_to_fit();

  v_weight.clear();
  v_weight.shrink_to_fit();
}
