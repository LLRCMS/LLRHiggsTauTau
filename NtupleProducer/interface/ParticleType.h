// enum for electron - muon -tau
//
#ifndef ParticleType_h
#define ParticleType_h

class ParticleType{
  public:
    enum particleType {
    MUON = 0,
    ELECTRON = 1,
    TAU =2
    };

    ParticleType() {}
    ~ParticleType(){}  
};

#endif
