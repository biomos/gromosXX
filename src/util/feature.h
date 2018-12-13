/**
 * @file feature.h
 * features and their cross checks
 */

#ifndef INCLUDED_FEATURE_H
#define INCLUDED_FEATURE_H

namespace util {
  /**
   * @class Feature
   * holds a feature like polarisation or parallelization
   */
  class Feature {
  public:
    /** 
     * Constructor
     */
    Feature() {}
    /**
     * Constructor with data
     */
    Feature(std::string &key, std::string &desc, bool active) :
        m_key(key), m_description(desc), m_active(active) {}
    /**
     * Constructor with data
     */
    Feature(const char* key, const char* desc, bool active) :
        m_key(key), m_description(desc), m_active(active) {}
    /**
     * copy constructor
     */
    Feature(const Feature &f) {
      key(f.key());
      description(f.description());
      is_active(f.is_active());
    }
    
    /**
     * const accessor to key
     */
    std::string key() const {
      return m_key;
    }    
    /**
     * accessor to key
     */
    std::string & key() {
      return m_key;
    }
    /**
     * accessor to key
     */
    void key(std::string const & value) {
      m_key = value;
    }
    /**
     * const accessor to description
     */
    std::string description() const {
      return m_description;
    }    
    /**
     * accessor to description
     */
    std::string & description() {
      return m_description;
    }
    /**
     * accessor to description
     */
    void description(std::string const & value) {
      m_description = value;
    }
    /**
     * const accessor to active
     */
    bool is_active() const {
      return m_active;
    }    
    /**
     * accessor to active
     */
    bool & is_active() {
      return m_active;
    }
    /**
     * accessor to active
     */
    void is_active(bool const & value) {
      m_active = value;
    }
    
  protected:
    /*
     * key that is used to unlock feature
     */
    std::string m_key;
    /*
     * describes what it does
     */
    std::string m_description;
    /*
     * determines whether the feature is active or not
     */
    bool m_active;
  };
  
  /**
   * @class FeatureChecker
   * checks features against each other
   */
  class FeatureChecker {
  public: 
    /**
     * @enum state
     */
    enum state {
      /**
       * feature locked
       */
      fc_locked,
      /**
       * feature unlocked
       */
      fc_unlocked,
      /**
       * feature unlocked with warning
       */
      fc_unlocked_warning
    };
    
    /**
     * constructor
     */
    FeatureChecker() {}
    /**
     * constructor with arguments
     */
    FeatureChecker(std::vector<Feature> &features) : m_features(features) {
      lock_all();
    }
   
    /**
     * accessor to the features
     */
    std::vector<Feature> features() {
      return m_features;
    }    

    /**
     * accessor to the features. This locks them all!
     */
    void features(std::vector<Feature> &value) {
      m_features = value;
      lock_all();
    }    
    
    /*
     * lock all features
     */
    void lock_all();
    /*
     * unlock two features against each other
     */
    void unlock(std::string key1, std::string key2, state s = fc_unlocked);
    /*
     * check features combinations
     */
    bool check();
    
  protected:
    /**
     * the features
     */
    std::vector<Feature> m_features;
    /*
     * the lock
     */
    std::map<std::string, std::map<std::string, state> > m_locked;
  };
}

#endif
