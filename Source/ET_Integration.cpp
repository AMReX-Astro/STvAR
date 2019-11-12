#include <ET_Integration.H>

namespace Variable
{
    amrex::Vector<std::string> names;

    void Initialize()
    {
        #include <ET_Integration_Variables.H>
    }
}