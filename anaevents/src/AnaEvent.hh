#ifndef __AnaEvent_hh__
#define __AnaEvent_hh__

#include <iostream>
#include <vector>

class AnaEvent
{
    public:
        AnaEvent(long int evid)
        {
            m_evid     = evid;
            m_flavor   = -1;
            m_topology = -1;
            m_reaction = -1;
            m_target   = -1;
            m_sample   = -1;
            m_signal   = false;
            m_sig_type = -1;
            m_true_evt = false;
            m_enu_true = -999.0;
            m_enu_reco = -999.0;
            m_trueD1   = -999.0;
            m_trueD2   = -999.0;
            m_recoD1   = -999.0;
            m_recoD2   = -999.0;
            m_q2_true  = -999.0;
            m_q2_reco  = -999.0;
            m_wght     = 1.0;
            m_wghtMC   = 1.0;
        }

        //Set/Get methods
        inline void SetTopology(const short val){ m_topology = val; }
        inline short GetTopology() const { return m_topology; }

        inline void SetReaction(const short val){ m_reaction = val; }
        inline short GetReaction() const { return m_reaction; }

        inline void SetTarget(const short val){ m_target = val; }
        inline short GetTarget() const { return m_target; }

        inline void SetSampleType(const short val){ m_sample = val; }
        inline short GetSampleType() const { return m_sample; }

        inline void SetSignalEvent(const bool flag = true){ m_signal = flag; }
        inline bool isSignalEvent() const { return m_signal; }

        inline void SetSignalType(const short val){ m_sig_type = val; }
        inline short GetSignalType() const { return m_sig_type; }

        inline void SetTrueEvent(const bool flag = true){ m_true_evt = flag; }
        inline bool isTrueEvent() const { return m_true_evt; }

        inline void SetFlavor(const short flavor){ m_flavor = flavor; }
        inline short GetFlavor() const { return m_flavor; }

        inline long int GetEvId() const { return m_evid; }

        inline void SetTrueEnu(double val) {m_enu_true = val;}
        inline double GetTrueEnu() const { return m_enu_true; }

        inline void SetRecoEnu(double val){ m_enu_reco = val; }
        inline double GetRecoEnu() const { return m_enu_reco; }

        inline void SetTrueD1(double val){ m_trueD1 = val; }
        inline double GetTrueD1() const { return m_trueD1; }

        inline void SetRecoD1(double val){ m_recoD1 = val; }
        inline double GetRecoD1() const { return m_recoD1; }

        inline void SetTrueD2(double val){ m_trueD2 = val; }
        inline double GetTrueD2() const { return m_trueD2; }

        inline void SetRecoD2(double val){ m_recoD2 = val; }
        inline double GetRecoD2() const { return m_recoD2; }

        inline void SetEvWght(double val){ m_wght  = val; }
        inline void SetEvWghtMC(double val){ m_wghtMC  = val; }
        inline void AddEvWght(double val){ m_wght *= val; }
        inline double GetEvWght() const { return m_wght; }
        inline double GetEvWghtMC() const { return m_wghtMC; }

        inline void ResetEvWght(){ m_wght = m_wghtMC; }

        inline void SetQ2Reco(double val){m_q2_reco = val;}
        inline double GetQ2Reco() const { return m_q2_reco; }

        inline void SetQ2True(double val){m_q2_true = val;}
        inline double GetQ2True() const { return m_q2_true; }

        void Print() const
        {
            std::cout << "Event ID    " << m_evid << std::endl
                      << "Topology    " << GetTopology() << std::endl
                      << "Reaction    " << GetReaction() << std::endl
                      << "Target      " << GetTarget() << std::endl
                      << "Flavor      " << GetFlavor() << std::endl
                      << "Sample      " << GetSampleType() << std::endl
                      << "Signal      " << GetSignalType() << std::endl
                      << "True energy " << GetTrueEnu() << std::endl
                      << "Reco energy " << GetRecoEnu() << std::endl
                      << "True D1     " << GetTrueD1() << std::endl
                      << "Reco D1     " << GetRecoD1() << std::endl
                      << "True D2     " << GetTrueD2() << std::endl
                      << "Reco D2     " << GetRecoD2() << std::endl
                      << "Weight      " << GetEvWght() << std::endl
                      << "Weight MC   " << GetEvWghtMC() << std::endl;
        }

        int GetEventVar(const std::string& var) const
        {
            if(var == "topology")
                return m_topology;
            else if(var == "reaction")
                return m_reaction;
            else if(var == "target")
                return m_target;
            else if(var == "flavor")
                return m_flavor;
            else if(var == "sample")
                return m_sample;
            else if(var == "signal")
                return m_sig_type;
            else
                return -1;
        }

        inline const std::vector<float>& GetRecoVar() const { return reco_var; }
        inline const std::vector<float>& GetTrueVar() const { return true_var; }
        inline void SetRecoVar(std::vector<float> vec) { reco_var = vec; }
        inline void SetTrueVar(std::vector<float> vec) { true_var = vec; }

    private:
        long int m_evid;   //unique event id
        short m_flavor;      //flavor of neutrino (numu, etc.)
        short m_topology;    //final state topology type
        short m_reaction;    //event interaction mode
        short m_sample;      //sample type (aka cutBranch)
        short m_sig_type;
        short m_target;      //target nuclei
        bool m_signal;     //flag if signal event
        bool m_true_evt;   //flag if true event
        double m_enu_true; //true nu energy
        double m_enu_reco; //recon nu energy
        double m_trueD1;   //true D1
        double m_trueD2;   //true D2
        double m_recoD1;   //reco D1
        double m_recoD2;   //reco D2
        double m_q2_true;
        double m_q2_reco;
        double m_wght;     //event weight
        double m_wghtMC;   //event weight from original MC

        //unsigned short nvar;
        std::vector<float> reco_var;
        std::vector<float> true_var;
};

#endif
