#ifndef OBSERVABLE_H
#define OBSERVABLE_H

class Observable {
public:
    virtual void Observe(const Configuration & conf, const ParamModel & param_mode ) = 0;
    virtual void SetZero() = 0;
    virtual void SetZeroAx() = 0;
};


#endif // OBSERVABLE_H
