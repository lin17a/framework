// -*- C++ -*-
// author: afiq anuar
// short: please refer to header for information

Framework::Tree::Tree(const std::string &filename, const std::string &treename, int compression /*= 505*/)
{
  file = std::make_unique<TFile>(filename.c_str(), "recreate", "", compression);
  file->cd();

  ptr = new TTree(treename.c_str(), "");
  ptr->SetAutoSave(0);
  ptr->SetImplicitMT(false);
}



template <typename Group, typename ...Attributes>
bool Framework::Tree::make_single_branches(Group &group, Attributes &&...attrs)
{
  auto iG = std::find_if(std::begin(v_branch), std::end(v_branch), 
                         [&name = group.name] (const auto &branch) {return name == std::get<0>(branch);});
  if (iG != std::end(v_branch))
    return false;

  auto iA = (group.has_attribute(attrs) and ...);
  if (!iA)
    throw std::invalid_argument( "ERROR: Tree::make_single_branches: some of the requested attributes are not within the group!!" );

  v_branch.emplace_back(group.name, nullptr, std::vector<TBranch *>(sizeof...(attrs), nullptr));

  std::array<std::string, sizeof...(attrs)> v_attr{ attrs... };
  for (int iB = 0; iB < v_attr.size(); ++iB) {
    std::visit(overload {
        // lolk apparently both "::" and "--" are invalid in branch names
        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB]] (std::vector<int> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), &vec[0], (name + "_" + attr + "/I").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB]] (std::vector<uint> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), &vec[0], (name + "_" + attr + "/i").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB]] (std::vector<float> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), &vec[0], (name + "_" + attr + "/F").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB]] (std::vector<double> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), &vec[0], (name + "_" + attr + "/D").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB]] (std::vector<long> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), &vec[0], (name + "_" + attr + "/L").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB]] (std::vector<ulong> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), &vec[0], (name + "_" + attr + "/l").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB]] (std::vector<boolean> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), &vec[0], (name + "_" + attr + "/O").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB]] (std::vector<char> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), &vec[0], (name + "_" + attr + "/B").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB]] (std::vector<unsigned char> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), &vec[0], (name + "_" + attr + "/b").c_str());},

        [] (auto &vec)
          { static_assert(always_false_v<decltype(vec)>, "ERROR: Tree::make_single_branches: unsupported type encountered!!"
                          " If it should have been supported, please add it and/or contact the developer."); }
      }, group.mref_to_attribute(v_attr[iB]));
  }

  std::get<1>(v_branch.back()) = std::function<void()>([&group] () { group.reorder(); });
  return true;
}



template <typename Group, typename ...Attributes>
bool Framework::Tree::make_array_branches(Group &group, Attributes &&...attrs)
{
  auto iG = std::find_if(std::begin(v_branch), std::end(v_branch), 
                         [&name = group.name] (const auto &branch) {return name == std::get<0>(branch);});
  if (iG != std::end(v_branch))
    return false;

  auto iA = (group.has_attribute(attrs) and ...);
  if (!iA)
    throw std::invalid_argument( "ERROR: Tree::make_array_branches: some of the requested attributes are not within the group!!" );

  v_branch.emplace_back(group.name, nullptr, std::vector<TBranch *>(sizeof...(attrs) + 1, nullptr));
  std::get<2>(v_branch.back())[0] = ptr->Branch(("n_" + group.name).c_str(), &group.mref_to_n_elements(), ("n_" + group.name + "/I").c_str());

  std::array<std::string, sizeof...(attrs)> v_attr{ attrs... };
  for (int iB = 0; iB < v_attr.size(); ++iB) {
    std::visit(overload {
        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB + 1]] (std::vector<int> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), vec.data(), (name + "_" + attr + "[n_"+ name + "]/I").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB + 1]] (std::vector<uint> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), vec.data(), (name + "_" + attr + "[n_"+ name + "]/i").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB + 1]] (std::vector<float> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), vec.data(), (name + "_" + attr + "[n_"+ name + "]/F").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB + 1]] (std::vector<double> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), vec.data(), (name + "_" + attr + "[n_"+ name + "]/D").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB + 1]] (std::vector<long> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), vec.data(), (name + "_" + attr + "[n_"+ name + "]/L").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB + 1]] (std::vector<ulong> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), vec.data(), (name + "_" + attr + "[n_"+ name + "]/l").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB + 1]] (std::vector<boolean> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), vec.data(), (name + "_" + attr + "[n_"+ name + "]/O").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB + 1]] (std::vector<char> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), vec.data(), (name + "_" + attr + "[n_"+ name + "]/B").c_str());},

        [&ptr = ptr, &attr = v_attr[iB], &name = group.name, &branch = std::get<2>(v_branch.back())[iB + 1]] (std::vector<unsigned char> &vec)
          {branch = ptr->Branch((name + "_" + attr).c_str(), vec.data(), (name + "_" + attr + "[n_"+ name + "]/b").c_str());},

        [] (auto &vec)
          { static_assert(always_false_v<decltype(vec)>, "ERROR: Tree::make_array_branches: unsupported type encountered!!"
                          " If it should have been supported, please add it and/or contact the developer."); }
      }, group.mref_to_attribute(v_attr[iB]));
  }

  int cap = group.ref_to_indices().capacity();
  std::get<1>(v_branch.back()) = std::function<void()>([&group, cap, v_attr, &branches = std::get<2>(v_branch.back())] () mutable {
      if (group.ref_to_indices().capacity() > cap) {
        for (int iB = 0; iB < v_attr.size(); ++iB)
          std::visit([&branch = branches[iB + 1]] (auto &vec) { branch->SetAddress(vec.data()); }, group.mref_to_attribute(v_attr[iB]));

        cap = group.ref_to_indices().capacity();
      }

      group.reorder();
    });
  return true;
}



void Framework::Tree::fill()
{
  for (auto &branch : v_branch)
    std::get<1>(branch)();

  ptr->Fill();
}



void Framework::Tree::save() const
{
  file->cd();
  ptr->Write();
}

