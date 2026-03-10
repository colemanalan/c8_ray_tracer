namespace c8_tracer
{

#define STOP_CLOSE_LENGTH 0.000005 // max distance for stopping criteria in Brent Loops

  /**
   * Utility function for Brent root finding searches. This generates and caches potential direction
   * and the function results
   */
  class CosineScanner
  {
  public:
    /**
     * @param cosMin minimum cosine value
     * @param cosMax maximum cosine value
     */
    CosineScanner(double cosMin, double cosMax)
        : cosMin_(cosMin), cosMax_(cosMax) {}

    CosineScanner() = delete;

    // check if the cosine value `x` has already been cached
    bool Has(double x) const
    {
      return index_.count(x) > 0;
    }

    // get the function result for cosine-value `x`
    LengthType Get(double x) const
    {
      return values_[index_.at(x)].second;
    }

    // store the function value `err` forcosine-value `x`
    void Cache(double x, LengthType err)
    {
      auto it = index_.find(x);
      if (it != index_.end())
      {
        values_[it->second].second = err;
        return;
      }
      index_[x] = values_.size();
      values_.emplace_back(x, err);
    }

    size_t CacheSize() const { return values_.size(); }

    // removes all invalid results (i.e. nan/inf from the cache)
    void RemoveIfInvalid()
    {
      // compact in-place
      size_t write = 0;
      for (size_t read = 0; read < values_.size(); ++read)
      {
        if (!std::isinf(values_[read].second))
        {
          values_[write] = values_[read];
          index_[values_[write].first] = write;
          write++;
        }
        else
        {
          index_.erase(values_[read].first);
        }
      }
      values_.resize(write);
    }

    // returns a list of sorted cosine values
    std::vector<double> SortedCosines() const
    {
      std::vector<double> xs;
      xs.reserve(values_.size());
      for (auto const &p : values_)
      {
        if (p.first >= cosMin_ && p.first <= cosMax_)
          xs.push_back(p.first);
      }
      std::sort(xs.begin(), xs.end());
      return xs;
    }

    // returns a vector of start-stop pairs corresponding to where the sign changes
    std::vector<std::pair<double, double>> FindSignChangeIntervals() const
    {
      auto xs = SortedCosines();
      std::vector<std::pair<double, double>> out;

      for (size_t i = 0; i + 1 < xs.size(); ++i)
      {
        double x0 = xs[i];
        double x1 = xs[i + 1];
        double f0 = Get(x0);
        double f1 = Get(x1);

        if (f0 * f1 <= 0.0 ||
            std::abs(f0) < STOP_CLOSE_LENGTH ||
            std::abs(f1) < STOP_CLOSE_LENGTH)
        {
          out.emplace_back(x0, x1);
        }
      }
      return out;
    }

    void UpdateBoundsFromBest(size_t nBest)
    {
      std::vector<std::pair<double, LengthType>> sorted = values_;
      std::sort(sorted.begin(), sorted.end(),
                [](auto const &a, auto const &b)
                {
                  return std::abs(a.second) < std::abs(b.second);
                });

      nBest = std::min(nBest, sorted.size());
      double lo = sorted[0].first;
      double hi = sorted[0].first;

      for (size_t i = 1; i < nBest; ++i)
      {
        lo = std::min(lo, sorted[i].first);
        hi = std::max(hi, sorted[i].first);
      }

      cosMin_ = std::max(-1.0, lo);
      cosMax_ = std::min(1.0, hi);
    }

    std::vector<double> GenerateGrid(int nRays) const
    {
      std::vector<double> xs(nRays);
      for (int i = 0; i < nRays; ++i)
      {
        double t = double(i) / double(nRays - 1);
        xs[i] = std::clamp(cosMin_ + t * (cosMax_ - cosMin_), -1.0, 1.0);
      }
      return xs;
    }

  private:
    double cosMin_;
    double cosMax_;
    std::vector<std::pair<double, LengthType>> values_;
    std::unordered_map<double, size_t> index_;
  };

} // namespace c8_tracer