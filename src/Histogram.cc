// -*- C++ -*-
// author: afiq anuar
// short: please refer to header for information

Framework::Histogram::Histogram()
{
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2(true);
}



template <typename Weighter>
void Framework::Histogram::set_weighter(Weighter weighter_)
{
  using Traits = function_traits<decltype(weighter_)>;
  static_assert(Traits::arity == 0, 
                "ERROR: Histogram::set_weighter: the number of arguments of the Histogram weighter must be zero."
                "Use lambda captures if some dependence on event information is needed.");
  static_assert(std::is_convertible_v<typename Traits::result_type, double>, 
                "ERROR: Histogram::set_weighter: the return type needs to be convertible to double!!");

  if (!weighter)
    weighter = std::function<double()>(weighter_);
}



template <typename Hist, typename Filler, typename ...Args>
bool Framework::Histogram::make_histogram(Filler filler, const std::string &name, Args &&...args)
{
  auto iH = std::find_if(std::begin(v_hist), std::end(v_hist), [&name] (const auto &hist) {return name == std::string(hist.first->GetName());});
  if (iH != std::end(v_hist))
    return false;

  auto f_fill = [this, index = v_hist.size(), filler] () {
    using Traits = function_traits<decltype(filler)>;

    filler((typename Traits::template bare_arg<0>) v_hist[index].first.get(), weight);
  };

  v_hist.emplace_back(std::make_unique<Hist>(name.c_str(), std::forward<Args>(args)...), std::function<void()>(f_fill));
  return true;
}



void Framework::Histogram::fill()
{
  weight = (weighter) ? weighter() : 1.;

  for (auto &hist : v_hist)
    hist.second();
}



void Framework::Histogram::unroll(bool add_under_overflow /*= true*/)
{
  for (int iH = 0; iH < v_hist.size(); ++iH) {
    std::string name = v_hist[iH].first->GetName();

    auto h2 = dynamic_cast<TH2 *>(v_hist[iH].first.get());
    auto h3 = dynamic_cast<TH3 *>(v_hist[iH].first.get());

    if (h2 != nullptr) {
      const int nx = h2->GetNbinsX(), ny = h2->GetNbinsY();
      v_hist.emplace_back(std::make_unique<TH1D>((name + "_unroll").c_str(), "", nx * ny, 0., nx * ny), std::function<void()>(nullptr));
      auto &hunr = v_hist.back().first;

      for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
          double content = h2->GetBinContent(ix + 1, iy + 1), squnc = h2->GetBinError(ix + 1, iy + 1) * h2->GetBinError(ix + 1, iy + 1);

          if (add_under_overflow) {
            if (ix == 0) {
              content += h2->GetBinContent(ix, iy + 1);
              squnc += h2->GetBinError(ix, iy + 1) * h2->GetBinError(ix, iy + 1);
            }
            if (iy == 0) {
              content += h2->GetBinContent(ix + 1, iy);
              squnc += h2->GetBinError(ix + 1, iy) * h2->GetBinError(ix + 1, iy);
            }
            if (ix == nx - 1) {
              content += h2->GetBinContent(ix + 2, iy + 1);
              squnc += h2->GetBinError(ix + 2, iy + 1) * h2->GetBinError(ix + 2, iy + 1);
            }
            if (iy == ny - 1) {
              content += h2->GetBinContent(ix + 1, iy + 2);
              squnc += h2->GetBinError(ix + 1, iy + 2) * h2->GetBinError(ix + 1, iy + 2);
            }
          }

          hunr->SetBinContent( (iy * nx) + ix + 1, content );
          hunr->SetBinError( (iy * nx) + ix + 1, std::sqrt(squnc) );
        }
      }
    }
    else if (h3 != nullptr) {
      const int nx = h3->GetNbinsX(), ny = h3->GetNbinsY(), nz = h3->GetNbinsZ();
      v_hist.emplace_back(std::make_unique<TH1D>((name + "_unroll").c_str(), "", nx * ny * nz, 0., nx * ny * nz), std::function<void()>(nullptr));
      auto &hunr = v_hist.back().first;

      for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
          for (int iz = 0; iz < nz; ++iz) {
            double content = h3->GetBinContent(ix + 1, iy + 1, iz + 1), 
              squnc = h3->GetBinError(ix + 1, iy + 1, iz + 1) * h3->GetBinError(ix + 1, iy + 1, iz + 1);

            if (add_under_overflow) {
              if (ix == 0) {
                content += h3->GetBinContent(ix, iy + 1, iz + 1);
                squnc += h3->GetBinError(ix, iy + 1, iz + 1) * h3->GetBinError(ix, iy + 1, iz + 1);
              }
              if (iy == 0) {
                content += h3->GetBinContent(ix + 1, iy, iz + 1);
                squnc += h3->GetBinError(ix + 1, iy, iz + 1) * h3->GetBinError(ix + 1, iy, iz + 1);
              }
              if (iz == 0) {
                content += h3->GetBinContent(ix + 1, iy + 1, iz);
                squnc += h3->GetBinError(ix + 1, iy + 1, iz) * h3->GetBinError(ix + 1, iy + 1, iz);
              }
              if (ix == nx - 1) {
                content += h3->GetBinContent(ix + 2, iy + 1, iz + 1);
                squnc += h3->GetBinError(ix + 2, iy + 1, iz + 1) * h3->GetBinError(ix + 2, iy + 1, iz + 1);
              }
              if (iy == ny - 1) {
                content += h3->GetBinContent(ix + 1, iy + 2, iz + 1);
                squnc += h3->GetBinError(ix + 1, iy + 2, iz + 1) * h3->GetBinError(ix + 1, iy + 2, iz + 1);
              }
              if (iz == nz - 1) {
                content += h3->GetBinContent(ix + 1, iy + 1, iz + 2);
                squnc += h3->GetBinError(ix + 1, iy + 1, iz + 2) * h3->GetBinError(ix + 1, iy + 1, iz + 2);
              }
            }

            hunr->SetBinContent( (iz * ny * nx) + (iy * nx) + ix + 1, content );
            hunr->SetBinError( (iz * ny * nx) + (iy * nx) + ix + 1, std::sqrt(squnc) );
          }
        }
      }
    }
  }
}



void Framework::Histogram::save_as(const std::string &name) const
{
  auto file = std::make_unique<TFile>(name.c_str(), "recreate");
  file->cd();

  for (auto &hist : v_hist)
    hist.first->Write();
}



const std::vector<Framework::Histogram::histfunc>& Framework::Histogram::histograms() const
{
  return v_hist;
}



template <typename ...Hists>
void Framework::save_all_as(const std::string &name, const Hists &...hists)
{
  using ref_to_vhf = std::reference_wrapper<const std::vector<Framework::Histogram::histfunc>>;
  std::array<ref_to_vhf, sizeof...(hists)> refs = { std::cref(hists.histograms())... };

  auto file = std::make_unique<TFile>(name.c_str(), "recreate");
  file->cd();

  for (const auto &ref : refs) {
    for (const auto &hist: ref.get())
      hist.first->Write();
  }
}
