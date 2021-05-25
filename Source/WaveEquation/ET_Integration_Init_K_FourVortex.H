#ifndef ET_INTEGRATION_INIT_K_H
#define ET_INTEGRATION_INIT_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_init(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab, 
        amrex::Real time, const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo(); 
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        
        amrex::Real Maxx = 20.0;
        
        amrex::Real v = 5;
        
        amrex::Real a = 3*2.33117;
        
        amrex::Real e = 0.3/std::sqrt(2);

        amrex::Real vboost = 0.5;

        amrex::Real gammaboost = std::pow(1 - std::pow(vboost, 2), -1.0/2.0);
        
        amrex::Real tm = -0.0001;

        amrex::Real tp = 0.0001;

        amrex::Real angle1_fidel = -M_PI/4;
        
        amrex::Real angle2_fidel = M_PI/4;
        
        amrex::Real angle1 = -angle2_fidel;
        
        amrex::Real angle2 = -angle1_fidel;

        amrex::Real ylprime = x1*std::cos(angle1)-x2*std::sin(angle1);

        amrex::Real yrprime = x1*std::cos(angle2)-x2*std::sin(angle2);
        
        amrex::Real xpl = gammaboost*(a - tp*vboost + x0);

        amrex::Real rpl = std::sqrt(std::pow(xpl, 2) + std::pow(ylprime, 2));

        amrex::Real xml = gammaboost*(a - tm*vboost + x0);

        amrex::Real rml = std::sqrt(std::pow(xml, 2) + std::pow(ylprime, 2));

        amrex::Real xpr = gammaboost*(-a + tp*vboost + x0);

        amrex::Real rpr = std::sqrt(std::pow(xpr, 2) + std::pow(yrprime, 2));

        amrex::Real xmr = gammaboost*(-a + tm*vboost + x0);

        amrex::Real rmr = std::sqrt(std::pow(xmr, 2) + std::pow(yrprime, 2));
        
        /*
        amrex::Real xpl = gammaboost*(a - tp*vboost + x0);

        amrex::Real rpl = std::sqrt(std::pow(xpl, 2) + std::pow(y, 2));

        amrex::Real xml = gammaboost*(a - tm*vboost + x0);

        amrex::Real rml = std::sqrt(std::pow(xml, 2) + std::pow(y, 2));

        amrex::Real xpr = gammaboost*(-a + tp*vboost + x0);

        amrex::Real rpr = std::sqrt(std::pow(xpr, 2) + std::pow(y, 2));

        amrex::Real xmr = gammaboost*(-a + tm*vboost + x0);

        amrex::Real rmr = std::sqrt(std::pow(xmr, 2) + std::pow(y, 2));
        */
        /*
        amrex::Real Fpl = 0;
        amrex::Real Fml = 0;
        amrex::Real Fpr = 0;
        amrex::Real Fmr = 0;
        amrex::Real Chipl = 0;
        amrex::Real Chiml = 0;
        amrex::Real Chipr = 0;
        amrex::Real Chimr = 0;
        */
        
        
        amrex::Real Fpl = v;
        amrex::Real Fml = v;
        amrex::Real Fpr = v;
        amrex::Real Fmr = v;
        amrex::Real Chipl = 0;
        amrex::Real Chiml = 0;
        amrex::Real Chipr = 0;
        amrex::Real Chimr = 0;
        amrex::Real apl = 0;
        amrex::Real aml = 0;
        amrex::Real apr = 0;
        amrex::Real amr = 0;

        if(rpl < Maxx) 
 { 

         Fpl = -1.0791825914970268*std::cos(2*std::acos(-1 + 2*rpl/Maxx)) + 0.80285889978116265*std::cos(3*std::acos(-1 + 2*rpl/Maxx)) - 0.54006781332116993*std::cos(4*std::acos(-1 + 2*rpl/Maxx)) + 0.33041822203594301*std::cos(5*std::acos(-1 + 2*rpl/Maxx)) - 0.18265025176158761*std::cos(6*std::acos(-1 + 2*rpl/Maxx)) + 0.088143335337428227*std::cos(7*std::acos(-1 + 2*rpl/Maxx)) - 0.032914921489218332*std::cos(8*std::acos(-1 + 2*rpl/Maxx)) + 0.0039265669028494682*std::cos(9*std::acos(-1 + 2*rpl/Maxx)) + 0.0088684564531565507*std::cos(10*std::acos(-1 + 2*rpl/Maxx)) - 0.012487917176595156*std::cos(11*std::acos(-1 + 2*rpl/Maxx)) + 0.011512829661162502*std::cos(12*std::acos(-1 + 2*rpl/Maxx)) - 0.008736996468339267*std::cos(13*std::acos(-1 + 2*rpl/Maxx)) + 0.0057139888532271449*std::cos(14*std::acos(-1 + 2*rpl/Maxx)) - 0.0031811020490695868*std::cos(15*std::acos(-1 + 2*rpl/Maxx)) + 0.0013766142424895406*std::cos(16*std::acos(-1 + 2*rpl/Maxx)) - 0.0002666432511588962*std::cos(17*std::acos(-1 + 2*rpl/Maxx)) - 0.00030193860802503469*std::cos(18*std::acos(-1 + 2*rpl/Maxx)) + 0.00050814209550134604*std::cos(19*std::acos(-1 + 2*rpl/Maxx)) - 0.00050768215750462975*std::cos(20*std::acos(-1 + 2*rpl/Maxx)) + 0.00041424714553496556*std::cos(21*std::acos(-1 + 2*rpl/Maxx)) - 0.00029887944196967109*std::cos(22*std::acos(-1 + 2*rpl/Maxx)) + 0.00019832092411475877*std::cos(23*std::acos(-1 + 2*rpl/Maxx)) - 0.00012591855763192938*std::cos(24*std::acos(-1 + 2*rpl/Maxx)) + 8.1570792144811111e-5*std::cos(25*std::acos(-1 + 2*rpl/Maxx)) - 5.9127074289927473e-5*std::cos(26*std::acos(-1 + 2*rpl/Maxx)) + 5.1061133080856913e-5*std::cos(27*std::acos(-1 + 2*rpl/Maxx)) - 5.0876555611187313e-5*std::cos(28*std::acos(-1 + 2*rpl/Maxx)) + 5.396483258066076e-5*std::cos(29*std::acos(-1 + 2*rpl/Maxx)) - 5.7571240465195683e-5*std::cos(30*std::acos(-1 + 2*rpl/Maxx)) + 6.0370580947335939e-5*std::cos(31*std::acos(-1 + 2*rpl/Maxx)) - 6.1961389745610052e-5*std::cos(32*std::acos(-1 + 2*rpl/Maxx)) + 6.2437584665536727e-5*std::cos(33*std::acos(-1 + 2*rpl/Maxx)) - 6.2089489739611259e-5*std::cos(34*std::acos(-1 + 2*rpl/Maxx)) + 6.1228909150320976e-5*std::cos(35*std::acos(-1 + 2*rpl/Maxx)) - 6.0107686167987229e-5*std::cos(36*std::acos(-1 + 2*rpl/Maxx)) + 5.8895261093941387e-5*std::cos(37*std::acos(-1 + 2*rpl/Maxx)) - 5.7687291804345286e-5*std::cos(38*std::acos(-1 + 2*rpl/Maxx)) + 5.6525168130025505e-5*std::cos(39*std::acos(-1 + 2*rpl/Maxx)) - 5.5417299469798231e-5*std::cos(40*std::acos(-1 + 2*rpl/Maxx)) + 5.4353656013819696e-5*std::cos(41*std::acos(-1 + 2*rpl/Maxx)) - 5.3320161792283494e-5*std::cos(42*std::acos(-1 + 2*rpl/Maxx)) + 5.2298314924745237e-5*std::cos(43*std::acos(-1 + 2*rpl/Maxx)) - 5.1285003905179871e-5*std::cos(44*std::acos(-1 + 2*rpl/Maxx)) + 5.0263916858819813e-5*std::cos(45*std::acos(-1 + 2*rpl/Maxx)) - 4.9283781304740315e-5*std::cos(46*std::acos(-1 + 2*rpl/Maxx)) + 4.8307125567192966e-5*std::cos(47*std::acos(-1 + 2*rpl/Maxx)) - 4.7592632780077415e-5*std::cos(48*std::acos(-1 + 2*rpl/Maxx)) + 2.3453717968781318e-5*std::cos(49*std::acos(-1 + 2*rpl/Maxx)) + 3.0117542335016725 + 2.5949803874590023*rpl/Maxx;


         Chipl = 1.1795075506129167*std::cos(2*std::acos(-1 + 2*rpl/Maxx)) - 0.53917011624401534*std::cos(3*std::acos(-1 + 2*rpl/Maxx)) + 0.16962696949652423*std::cos(4*std::acos(-1 + 2*rpl/Maxx)) - 0.0037342126895790538*std::cos(5*std::acos(-1 + 2*rpl/Maxx)) - 0.053759834803502007*std::cos(6*std::acos(-1 + 2*rpl/Maxx)) + 0.062073415559850888*std::cos(7*std::acos(-1 + 2*rpl/Maxx)) - 0.051729276271027136*std::cos(8*std::acos(-1 + 2*rpl/Maxx)) + 0.036745573212263516*std::cos(9*std::acos(-1 + 2*rpl/Maxx)) - 0.023006732885378415*std::cos(10*std::acos(-1 + 2*rpl/Maxx)) + 0.012502022934543887*std::cos(11*std::acos(-1 + 2*rpl/Maxx)) - 0.0054260815940578196*std::cos(12*std::acos(-1 + 2*rpl/Maxx)) + 0.0012235286986818584*std::cos(13*std::acos(-1 + 2*rpl/Maxx)) + 0.00088379036280092978*std::cos(14*std::acos(-1 + 2*rpl/Maxx)) - 0.0016379461520883594*std::cos(15*std::acos(-1 + 2*rpl/Maxx)) + 0.0016334052552107966*std::cos(16*std::acos(-1 + 2*rpl/Maxx)) - 0.0012885284225685094*std::cos(17*std::acos(-1 + 2*rpl/Maxx)) + 0.00086167705025049965*std::cos(18*std::acos(-1 + 2*rpl/Maxx)) - 0.00048717787299036548*std::cos(19*std::acos(-1 + 2*rpl/Maxx)) + 0.00021535150786820636*std::cos(20*std::acos(-1 + 2*rpl/Maxx)) - 4.7575985736512345e-5*std::cos(21*std::acos(-1 + 2*rpl/Maxx)) - 3.7538066526578624e-5*std::cos(22*std::acos(-1 + 2*rpl/Maxx)) + 6.7304605045776601e-5*std::cos(23*std::acos(-1 + 2*rpl/Maxx)) - 6.5929427728895202e-5*std::cos(24*std::acos(-1 + 2*rpl/Maxx)) + 5.1011723657978163e-5*std::cos(25*std::acos(-1 + 2*rpl/Maxx)) - 3.3360399039120973e-5*std::cos(26*std::acos(-1 + 2*rpl/Maxx)) + 1.8380576200193871e-5*std::cos(27*std::acos(-1 + 2*rpl/Maxx)) - 7.878931024914073e-6*std::cos(28*std::acos(-1 + 2*rpl/Maxx)) + 1.6644334001256316e-6*std::cos(29*std::acos(-1 + 2*rpl/Maxx)) + 1.2960757659209459e-6*std::cos(30*std::acos(-1 + 2*rpl/Maxx)) - 2.178037631388182e-6*std::cos(31*std::acos(-1 + 2*rpl/Maxx)) + 1.9581887788229492e-6*std::cos(32*std::acos(-1 + 2*rpl/Maxx)) - 1.3087487030288237e-6*std::cos(33*std::acos(-1 + 2*rpl/Maxx)) + 6.1900318628474647e-7*std::cos(34*std::acos(-1 + 2*rpl/Maxx)) - 6.8004820184653329e-8*std::cos(35*std::acos(-1 + 2*rpl/Maxx)) - 2.9665542454353737e-7*std::cos(36*std::acos(-1 + 2*rpl/Maxx)) + 4.9510036131430351e-7*std::cos(37*std::acos(-1 + 2*rpl/Maxx)) - 5.7227885640485122e-7*std::cos(38*std::acos(-1 + 2*rpl/Maxx)) + 5.7382153603081516e-7*std::cos(39*std::acos(-1 + 2*rpl/Maxx)) - 5.3555876621860833e-7*std::cos(40*std::acos(-1 + 2*rpl/Maxx)) + 4.8152687055446928e-7*std::cos(41*std::acos(-1 + 2*rpl/Maxx)) - 4.2637502776406691e-7*std::cos(42*std::acos(-1 + 2*rpl/Maxx)) + 3.795025697779657e-7*std::cos(43*std::acos(-1 + 2*rpl/Maxx)) - 3.4898785130852968e-7*std::cos(44*std::acos(-1 + 2*rpl/Maxx)) + 3.4383928358246674e-7*std::cos(45*std::acos(-1 + 2*rpl/Maxx)) - 3.7153750839260132e-7*std::cos(46*std::acos(-1 + 2*rpl/Maxx)) + 4.2910669084185182e-7*std::cos(47*std::acos(-1 + 2*rpl/Maxx)) - 4.9057561744160653e-7*std::cos(48*std::acos(-1 + 2*rpl/Maxx)) - 0.00024358893649518274*std::cos(49*std::acos(-1 + 2*rpl/Maxx)) + 3.14211374001962 - 3.9270934315856834*rpl/Maxx;


         apl = 0.26168769006825793*std::cos(2*std::acos(-1 + 2*rpl/Maxx)) - 0.18455890149887749*std::cos(3*std::acos(-1 + 2*rpl/Maxx)) + 0.10644691856087608*std::cos(4*std::acos(-1 + 2*rpl/Maxx)) - 0.042626922046416882*std::cos(5*std::acos(-1 + 2*rpl/Maxx)) + 0.00056248663002799262*std::cos(6*std::acos(-1 + 2*rpl/Maxx)) + 0.02014476610478377*std::cos(7*std::acos(-1 + 2*rpl/Maxx)) - 0.024957649545127161*std::cos(8*std::acos(-1 + 2*rpl/Maxx)) + 0.020849460083631154*std::cos(9*std::acos(-1 + 2*rpl/Maxx)) - 0.013666027854196879*std::cos(10*std::acos(-1 + 2*rpl/Maxx)) + 0.0070116747178595633*std::cos(11*std::acos(-1 + 2*rpl/Maxx)) - 0.0023727188691233798*std::cos(12*std::acos(-1 + 2*rpl/Maxx)) - 0.00016454816727440303*std::cos(13*std::acos(-1 + 2*rpl/Maxx)) + 0.0011512773593783492*std::cos(14*std::acos(-1 + 2*rpl/Maxx)) - 0.0012464661704612013*std::cos(15*std::acos(-1 + 2*rpl/Maxx)) + 0.00096190009053586376*std::cos(16*std::acos(-1 + 2*rpl/Maxx)) - 0.00060353310259851785*std::cos(17*std::acos(-1 + 2*rpl/Maxx)) + 0.00031014612973111227*std::cos(18*std::acos(-1 + 2*rpl/Maxx)) - 0.00011760298419894008*std::cos(19*std::acos(-1 + 2*rpl/Maxx)) + 1.2972426905246089e-5*std::cos(20*std::acos(-1 + 2*rpl/Maxx)) + 3.1460299297518064e-5*std::cos(21*std::acos(-1 + 2*rpl/Maxx)) - 4.1565148757710908e-5*std::cos(22*std::acos(-1 + 2*rpl/Maxx)) + 3.5796814181932425e-5*std::cos(23*std::acos(-1 + 2*rpl/Maxx)) - 2.5140448974753429e-5*std::cos(24*std::acos(-1 + 2*rpl/Maxx)) + 1.5034014273756846e-5*std::cos(25*std::acos(-1 + 2*rpl/Maxx)) - 7.481343748544104e-6*std::cos(26*std::acos(-1 + 2*rpl/Maxx)) + 2.6858337195837846e-6*std::cos(27*std::acos(-1 + 2*rpl/Maxx)) - 9.9747609096360928e-8*std::cos(28*std::acos(-1 + 2*rpl/Maxx)) - 9.9720849232281744e-7*std::cos(29*std::acos(-1 + 2*rpl/Maxx)) + 1.2325462535254711e-6*std::cos(30*std::acos(-1 + 2*rpl/Maxx)) - 1.0565043162282929e-6*std::cos(31*std::acos(-1 + 2*rpl/Maxx)) + 7.4718662447722397e-7*std::cos(32*std::acos(-1 + 2*rpl/Maxx)) - 4.503170181331782e-7*std::cos(33*std::acos(-1 + 2*rpl/Maxx)) + 2.2444015547850691e-7*std::cos(34*std::acos(-1 + 2*rpl/Maxx)) - 7.8667430662813549e-8*std::cos(35*std::acos(-1 + 2*rpl/Maxx)) - 7.5422423704999825e-10*std::cos(36*std::acos(-1 + 2*rpl/Maxx)) + 3.4257847261451717e-8*std::cos(37*std::acos(-1 + 2*rpl/Maxx)) - 4.0643714033846369e-8*std::cos(38*std::acos(-1 + 2*rpl/Maxx)) + 3.3910631625957225e-8*std::cos(39*std::acos(-1 + 2*rpl/Maxx)) - 2.2969828893536417e-8*std::cos(40*std::acos(-1 + 2*rpl/Maxx)) + 1.26468956304575e-8*std::cos(41*std::acos(-1 + 2*rpl/Maxx)) - 4.9989372907787308e-9*std::cos(42*std::acos(-1 + 2*rpl/Maxx)) + 4.542935982275102e-10*std::cos(43*std::acos(-1 + 2*rpl/Maxx)) + 1.3930406250537502e-9*std::cos(44*std::acos(-1 + 2*rpl/Maxx)) - 1.3113806317392608e-9*std::cos(45*std::acos(-1 + 2*rpl/Maxx)) + 1.7380610306497266e-10*std::cos(46*std::acos(-1 + 2*rpl/Maxx)) + 1.2009518158304569e-9*std::cos(47*std::acos(-1 + 2*rpl/Maxx)) - 2.1730906049510922e-9*std::cos(48*std::acos(-1 + 2*rpl/Maxx)) - 5.8801492099470041e-5*std::cos(49*std::acos(-1 + 2*rpl/Maxx)) + 0.48870556108500007 - 0.63754080718652051*rpl/Maxx;


}
if(rml < Maxx) 
 { 

         Fml = -1.0791825914970268*std::cos(2*std::acos(-1 + 2*rml/Maxx)) + 0.80285889978116265*std::cos(3*std::acos(-1 + 2*rml/Maxx)) - 0.54006781332116993*std::cos(4*std::acos(-1 + 2*rml/Maxx)) + 0.33041822203594301*std::cos(5*std::acos(-1 + 2*rml/Maxx)) - 0.18265025176158761*std::cos(6*std::acos(-1 + 2*rml/Maxx)) + 0.088143335337428227*std::cos(7*std::acos(-1 + 2*rml/Maxx)) - 0.032914921489218332*std::cos(8*std::acos(-1 + 2*rml/Maxx)) + 0.0039265669028494682*std::cos(9*std::acos(-1 + 2*rml/Maxx)) + 0.0088684564531565507*std::cos(10*std::acos(-1 + 2*rml/Maxx)) - 0.012487917176595156*std::cos(11*std::acos(-1 + 2*rml/Maxx)) + 0.011512829661162502*std::cos(12*std::acos(-1 + 2*rml/Maxx)) - 0.008736996468339267*std::cos(13*std::acos(-1 + 2*rml/Maxx)) + 0.0057139888532271449*std::cos(14*std::acos(-1 + 2*rml/Maxx)) - 0.0031811020490695868*std::cos(15*std::acos(-1 + 2*rml/Maxx)) + 0.0013766142424895406*std::cos(16*std::acos(-1 + 2*rml/Maxx)) - 0.0002666432511588962*std::cos(17*std::acos(-1 + 2*rml/Maxx)) - 0.00030193860802503469*std::cos(18*std::acos(-1 + 2*rml/Maxx)) + 0.00050814209550134604*std::cos(19*std::acos(-1 + 2*rml/Maxx)) - 0.00050768215750462975*std::cos(20*std::acos(-1 + 2*rml/Maxx)) + 0.00041424714553496556*std::cos(21*std::acos(-1 + 2*rml/Maxx)) - 0.00029887944196967109*std::cos(22*std::acos(-1 + 2*rml/Maxx)) + 0.00019832092411475877*std::cos(23*std::acos(-1 + 2*rml/Maxx)) - 0.00012591855763192938*std::cos(24*std::acos(-1 + 2*rml/Maxx)) + 8.1570792144811111e-5*std::cos(25*std::acos(-1 + 2*rml/Maxx)) - 5.9127074289927473e-5*std::cos(26*std::acos(-1 + 2*rml/Maxx)) + 5.1061133080856913e-5*std::cos(27*std::acos(-1 + 2*rml/Maxx)) - 5.0876555611187313e-5*std::cos(28*std::acos(-1 + 2*rml/Maxx)) + 5.396483258066076e-5*std::cos(29*std::acos(-1 + 2*rml/Maxx)) - 5.7571240465195683e-5*std::cos(30*std::acos(-1 + 2*rml/Maxx)) + 6.0370580947335939e-5*std::cos(31*std::acos(-1 + 2*rml/Maxx)) - 6.1961389745610052e-5*std::cos(32*std::acos(-1 + 2*rml/Maxx)) + 6.2437584665536727e-5*std::cos(33*std::acos(-1 + 2*rml/Maxx)) - 6.2089489739611259e-5*std::cos(34*std::acos(-1 + 2*rml/Maxx)) + 6.1228909150320976e-5*std::cos(35*std::acos(-1 + 2*rml/Maxx)) - 6.0107686167987229e-5*std::cos(36*std::acos(-1 + 2*rml/Maxx)) + 5.8895261093941387e-5*std::cos(37*std::acos(-1 + 2*rml/Maxx)) - 5.7687291804345286e-5*std::cos(38*std::acos(-1 + 2*rml/Maxx)) + 5.6525168130025505e-5*std::cos(39*std::acos(-1 + 2*rml/Maxx)) - 5.5417299469798231e-5*std::cos(40*std::acos(-1 + 2*rml/Maxx)) + 5.4353656013819696e-5*std::cos(41*std::acos(-1 + 2*rml/Maxx)) - 5.3320161792283494e-5*std::cos(42*std::acos(-1 + 2*rml/Maxx)) + 5.2298314924745237e-5*std::cos(43*std::acos(-1 + 2*rml/Maxx)) - 5.1285003905179871e-5*std::cos(44*std::acos(-1 + 2*rml/Maxx)) + 5.0263916858819813e-5*std::cos(45*std::acos(-1 + 2*rml/Maxx)) - 4.9283781304740315e-5*std::cos(46*std::acos(-1 + 2*rml/Maxx)) + 4.8307125567192966e-5*std::cos(47*std::acos(-1 + 2*rml/Maxx)) - 4.7592632780077415e-5*std::cos(48*std::acos(-1 + 2*rml/Maxx)) + 2.3453717968781318e-5*std::cos(49*std::acos(-1 + 2*rml/Maxx)) + 3.0117542335016725 + 2.5949803874590023*rml/Maxx;


         Chiml = 1.1795075506129167*std::cos(2*std::acos(-1 + 2*rml/Maxx)) - 0.53917011624401534*std::cos(3*std::acos(-1 + 2*rml/Maxx)) + 0.16962696949652423*std::cos(4*std::acos(-1 + 2*rml/Maxx)) - 0.0037342126895790538*std::cos(5*std::acos(-1 + 2*rml/Maxx)) - 0.053759834803502007*std::cos(6*std::acos(-1 + 2*rml/Maxx)) + 0.062073415559850888*std::cos(7*std::acos(-1 + 2*rml/Maxx)) - 0.051729276271027136*std::cos(8*std::acos(-1 + 2*rml/Maxx)) + 0.036745573212263516*std::cos(9*std::acos(-1 + 2*rml/Maxx)) - 0.023006732885378415*std::cos(10*std::acos(-1 + 2*rml/Maxx)) + 0.012502022934543887*std::cos(11*std::acos(-1 + 2*rml/Maxx)) - 0.0054260815940578196*std::cos(12*std::acos(-1 + 2*rml/Maxx)) + 0.0012235286986818584*std::cos(13*std::acos(-1 + 2*rml/Maxx)) + 0.00088379036280092978*std::cos(14*std::acos(-1 + 2*rml/Maxx)) - 0.0016379461520883594*std::cos(15*std::acos(-1 + 2*rml/Maxx)) + 0.0016334052552107966*std::cos(16*std::acos(-1 + 2*rml/Maxx)) - 0.0012885284225685094*std::cos(17*std::acos(-1 + 2*rml/Maxx)) + 0.00086167705025049965*std::cos(18*std::acos(-1 + 2*rml/Maxx)) - 0.00048717787299036548*std::cos(19*std::acos(-1 + 2*rml/Maxx)) + 0.00021535150786820636*std::cos(20*std::acos(-1 + 2*rml/Maxx)) - 4.7575985736512345e-5*std::cos(21*std::acos(-1 + 2*rml/Maxx)) - 3.7538066526578624e-5*std::cos(22*std::acos(-1 + 2*rml/Maxx)) + 6.7304605045776601e-5*std::cos(23*std::acos(-1 + 2*rml/Maxx)) - 6.5929427728895202e-5*std::cos(24*std::acos(-1 + 2*rml/Maxx)) + 5.1011723657978163e-5*std::cos(25*std::acos(-1 + 2*rml/Maxx)) - 3.3360399039120973e-5*std::cos(26*std::acos(-1 + 2*rml/Maxx)) + 1.8380576200193871e-5*std::cos(27*std::acos(-1 + 2*rml/Maxx)) - 7.878931024914073e-6*std::cos(28*std::acos(-1 + 2*rml/Maxx)) + 1.6644334001256316e-6*std::cos(29*std::acos(-1 + 2*rml/Maxx)) + 1.2960757659209459e-6*std::cos(30*std::acos(-1 + 2*rml/Maxx)) - 2.178037631388182e-6*std::cos(31*std::acos(-1 + 2*rml/Maxx)) + 1.9581887788229492e-6*std::cos(32*std::acos(-1 + 2*rml/Maxx)) - 1.3087487030288237e-6*std::cos(33*std::acos(-1 + 2*rml/Maxx)) + 6.1900318628474647e-7*std::cos(34*std::acos(-1 + 2*rml/Maxx)) - 6.8004820184653329e-8*std::cos(35*std::acos(-1 + 2*rml/Maxx)) - 2.9665542454353737e-7*std::cos(36*std::acos(-1 + 2*rml/Maxx)) + 4.9510036131430351e-7*std::cos(37*std::acos(-1 + 2*rml/Maxx)) - 5.7227885640485122e-7*std::cos(38*std::acos(-1 + 2*rml/Maxx)) + 5.7382153603081516e-7*std::cos(39*std::acos(-1 + 2*rml/Maxx)) - 5.3555876621860833e-7*std::cos(40*std::acos(-1 + 2*rml/Maxx)) + 4.8152687055446928e-7*std::cos(41*std::acos(-1 + 2*rml/Maxx)) - 4.2637502776406691e-7*std::cos(42*std::acos(-1 + 2*rml/Maxx)) + 3.795025697779657e-7*std::cos(43*std::acos(-1 + 2*rml/Maxx)) - 3.4898785130852968e-7*std::cos(44*std::acos(-1 + 2*rml/Maxx)) + 3.4383928358246674e-7*std::cos(45*std::acos(-1 + 2*rml/Maxx)) - 3.7153750839260132e-7*std::cos(46*std::acos(-1 + 2*rml/Maxx)) + 4.2910669084185182e-7*std::cos(47*std::acos(-1 + 2*rml/Maxx)) - 4.9057561744160653e-7*std::cos(48*std::acos(-1 + 2*rml/Maxx)) - 0.00024358893649518274*std::cos(49*std::acos(-1 + 2*rml/Maxx)) + 3.14211374001962 - 3.9270934315856834*rml/Maxx;


         aml = 0.26168769006825793*std::cos(2*std::acos(-1 + 2*rml/Maxx)) - 0.18455890149887749*std::cos(3*std::acos(-1 + 2*rml/Maxx)) + 0.10644691856087608*std::cos(4*std::acos(-1 + 2*rml/Maxx)) - 0.042626922046416882*std::cos(5*std::acos(-1 + 2*rml/Maxx)) + 0.00056248663002799262*std::cos(6*std::acos(-1 + 2*rml/Maxx)) + 0.02014476610478377*std::cos(7*std::acos(-1 + 2*rml/Maxx)) - 0.024957649545127161*std::cos(8*std::acos(-1 + 2*rml/Maxx)) + 0.020849460083631154*std::cos(9*std::acos(-1 + 2*rml/Maxx)) - 0.013666027854196879*std::cos(10*std::acos(-1 + 2*rml/Maxx)) + 0.0070116747178595633*std::cos(11*std::acos(-1 + 2*rml/Maxx)) - 0.0023727188691233798*std::cos(12*std::acos(-1 + 2*rml/Maxx)) - 0.00016454816727440303*std::cos(13*std::acos(-1 + 2*rml/Maxx)) + 0.0011512773593783492*std::cos(14*std::acos(-1 + 2*rml/Maxx)) - 0.0012464661704612013*std::cos(15*std::acos(-1 + 2*rml/Maxx)) + 0.00096190009053586376*std::cos(16*std::acos(-1 + 2*rml/Maxx)) - 0.00060353310259851785*std::cos(17*std::acos(-1 + 2*rml/Maxx)) + 0.00031014612973111227*std::cos(18*std::acos(-1 + 2*rml/Maxx)) - 0.00011760298419894008*std::cos(19*std::acos(-1 + 2*rml/Maxx)) + 1.2972426905246089e-5*std::cos(20*std::acos(-1 + 2*rml/Maxx)) + 3.1460299297518064e-5*std::cos(21*std::acos(-1 + 2*rml/Maxx)) - 4.1565148757710908e-5*std::cos(22*std::acos(-1 + 2*rml/Maxx)) + 3.5796814181932425e-5*std::cos(23*std::acos(-1 + 2*rml/Maxx)) - 2.5140448974753429e-5*std::cos(24*std::acos(-1 + 2*rml/Maxx)) + 1.5034014273756846e-5*std::cos(25*std::acos(-1 + 2*rml/Maxx)) - 7.481343748544104e-6*std::cos(26*std::acos(-1 + 2*rml/Maxx)) + 2.6858337195837846e-6*std::cos(27*std::acos(-1 + 2*rml/Maxx)) - 9.9747609096360928e-8*std::cos(28*std::acos(-1 + 2*rml/Maxx)) - 9.9720849232281744e-7*std::cos(29*std::acos(-1 + 2*rml/Maxx)) + 1.2325462535254711e-6*std::cos(30*std::acos(-1 + 2*rml/Maxx)) - 1.0565043162282929e-6*std::cos(31*std::acos(-1 + 2*rml/Maxx)) + 7.4718662447722397e-7*std::cos(32*std::acos(-1 + 2*rml/Maxx)) - 4.503170181331782e-7*std::cos(33*std::acos(-1 + 2*rml/Maxx)) + 2.2444015547850691e-7*std::cos(34*std::acos(-1 + 2*rml/Maxx)) - 7.8667430662813549e-8*std::cos(35*std::acos(-1 + 2*rml/Maxx)) - 7.5422423704999825e-10*std::cos(36*std::acos(-1 + 2*rml/Maxx)) + 3.4257847261451717e-8*std::cos(37*std::acos(-1 + 2*rml/Maxx)) - 4.0643714033846369e-8*std::cos(38*std::acos(-1 + 2*rml/Maxx)) + 3.3910631625957225e-8*std::cos(39*std::acos(-1 + 2*rml/Maxx)) - 2.2969828893536417e-8*std::cos(40*std::acos(-1 + 2*rml/Maxx)) + 1.26468956304575e-8*std::cos(41*std::acos(-1 + 2*rml/Maxx)) - 4.9989372907787308e-9*std::cos(42*std::acos(-1 + 2*rml/Maxx)) + 4.542935982275102e-10*std::cos(43*std::acos(-1 + 2*rml/Maxx)) + 1.3930406250537502e-9*std::cos(44*std::acos(-1 + 2*rml/Maxx)) - 1.3113806317392608e-9*std::cos(45*std::acos(-1 + 2*rml/Maxx)) + 1.7380610306497266e-10*std::cos(46*std::acos(-1 + 2*rml/Maxx)) + 1.2009518158304569e-9*std::cos(47*std::acos(-1 + 2*rml/Maxx)) - 2.1730906049510922e-9*std::cos(48*std::acos(-1 + 2*rml/Maxx)) - 5.8801492099470041e-5*std::cos(49*std::acos(-1 + 2*rml/Maxx)) + 0.48870556108500007 - 0.63754080718652051*rml/Maxx;


}
if(rpr < Maxx) 
 { 

         Fpr = -1.0791825914970268*std::cos(2*std::acos(-1 + 2*rpr/Maxx)) + 0.80285889978116265*std::cos(3*std::acos(-1 + 2*rpr/Maxx)) - 0.54006781332116993*std::cos(4*std::acos(-1 + 2*rpr/Maxx)) + 0.33041822203594301*std::cos(5*std::acos(-1 + 2*rpr/Maxx)) - 0.18265025176158761*std::cos(6*std::acos(-1 + 2*rpr/Maxx)) + 0.088143335337428227*std::cos(7*std::acos(-1 + 2*rpr/Maxx)) - 0.032914921489218332*std::cos(8*std::acos(-1 + 2*rpr/Maxx)) + 0.0039265669028494682*std::cos(9*std::acos(-1 + 2*rpr/Maxx)) + 0.0088684564531565507*std::cos(10*std::acos(-1 + 2*rpr/Maxx)) - 0.012487917176595156*std::cos(11*std::acos(-1 + 2*rpr/Maxx)) + 0.011512829661162502*std::cos(12*std::acos(-1 + 2*rpr/Maxx)) - 0.008736996468339267*std::cos(13*std::acos(-1 + 2*rpr/Maxx)) + 0.0057139888532271449*std::cos(14*std::acos(-1 + 2*rpr/Maxx)) - 0.0031811020490695868*std::cos(15*std::acos(-1 + 2*rpr/Maxx)) + 0.0013766142424895406*std::cos(16*std::acos(-1 + 2*rpr/Maxx)) - 0.0002666432511588962*std::cos(17*std::acos(-1 + 2*rpr/Maxx)) - 0.00030193860802503469*std::cos(18*std::acos(-1 + 2*rpr/Maxx)) + 0.00050814209550134604*std::cos(19*std::acos(-1 + 2*rpr/Maxx)) - 0.00050768215750462975*std::cos(20*std::acos(-1 + 2*rpr/Maxx)) + 0.00041424714553496556*std::cos(21*std::acos(-1 + 2*rpr/Maxx)) - 0.00029887944196967109*std::cos(22*std::acos(-1 + 2*rpr/Maxx)) + 0.00019832092411475877*std::cos(23*std::acos(-1 + 2*rpr/Maxx)) - 0.00012591855763192938*std::cos(24*std::acos(-1 + 2*rpr/Maxx)) + 8.1570792144811111e-5*std::cos(25*std::acos(-1 + 2*rpr/Maxx)) - 5.9127074289927473e-5*std::cos(26*std::acos(-1 + 2*rpr/Maxx)) + 5.1061133080856913e-5*std::cos(27*std::acos(-1 + 2*rpr/Maxx)) - 5.0876555611187313e-5*std::cos(28*std::acos(-1 + 2*rpr/Maxx)) + 5.396483258066076e-5*std::cos(29*std::acos(-1 + 2*rpr/Maxx)) - 5.7571240465195683e-5*std::cos(30*std::acos(-1 + 2*rpr/Maxx)) + 6.0370580947335939e-5*std::cos(31*std::acos(-1 + 2*rpr/Maxx)) - 6.1961389745610052e-5*std::cos(32*std::acos(-1 + 2*rpr/Maxx)) + 6.2437584665536727e-5*std::cos(33*std::acos(-1 + 2*rpr/Maxx)) - 6.2089489739611259e-5*std::cos(34*std::acos(-1 + 2*rpr/Maxx)) + 6.1228909150320976e-5*std::cos(35*std::acos(-1 + 2*rpr/Maxx)) - 6.0107686167987229e-5*std::cos(36*std::acos(-1 + 2*rpr/Maxx)) + 5.8895261093941387e-5*std::cos(37*std::acos(-1 + 2*rpr/Maxx)) - 5.7687291804345286e-5*std::cos(38*std::acos(-1 + 2*rpr/Maxx)) + 5.6525168130025505e-5*std::cos(39*std::acos(-1 + 2*rpr/Maxx)) - 5.5417299469798231e-5*std::cos(40*std::acos(-1 + 2*rpr/Maxx)) + 5.4353656013819696e-5*std::cos(41*std::acos(-1 + 2*rpr/Maxx)) - 5.3320161792283494e-5*std::cos(42*std::acos(-1 + 2*rpr/Maxx)) + 5.2298314924745237e-5*std::cos(43*std::acos(-1 + 2*rpr/Maxx)) - 5.1285003905179871e-5*std::cos(44*std::acos(-1 + 2*rpr/Maxx)) + 5.0263916858819813e-5*std::cos(45*std::acos(-1 + 2*rpr/Maxx)) - 4.9283781304740315e-5*std::cos(46*std::acos(-1 + 2*rpr/Maxx)) + 4.8307125567192966e-5*std::cos(47*std::acos(-1 + 2*rpr/Maxx)) - 4.7592632780077415e-5*std::cos(48*std::acos(-1 + 2*rpr/Maxx)) + 2.3453717968781318e-5*std::cos(49*std::acos(-1 + 2*rpr/Maxx)) + 3.0117542335016725 + 2.5949803874590023*rpr/Maxx;


         Chipr = 1.1795075506129167*std::cos(2*std::acos(-1 + 2*rpr/Maxx)) - 0.53917011624401534*std::cos(3*std::acos(-1 + 2*rpr/Maxx)) + 0.16962696949652423*std::cos(4*std::acos(-1 + 2*rpr/Maxx)) - 0.0037342126895790538*std::cos(5*std::acos(-1 + 2*rpr/Maxx)) - 0.053759834803502007*std::cos(6*std::acos(-1 + 2*rpr/Maxx)) + 0.062073415559850888*std::cos(7*std::acos(-1 + 2*rpr/Maxx)) - 0.051729276271027136*std::cos(8*std::acos(-1 + 2*rpr/Maxx)) + 0.036745573212263516*std::cos(9*std::acos(-1 + 2*rpr/Maxx)) - 0.023006732885378415*std::cos(10*std::acos(-1 + 2*rpr/Maxx)) + 0.012502022934543887*std::cos(11*std::acos(-1 + 2*rpr/Maxx)) - 0.0054260815940578196*std::cos(12*std::acos(-1 + 2*rpr/Maxx)) + 0.0012235286986818584*std::cos(13*std::acos(-1 + 2*rpr/Maxx)) + 0.00088379036280092978*std::cos(14*std::acos(-1 + 2*rpr/Maxx)) - 0.0016379461520883594*std::cos(15*std::acos(-1 + 2*rpr/Maxx)) + 0.0016334052552107966*std::cos(16*std::acos(-1 + 2*rpr/Maxx)) - 0.0012885284225685094*std::cos(17*std::acos(-1 + 2*rpr/Maxx)) + 0.00086167705025049965*std::cos(18*std::acos(-1 + 2*rpr/Maxx)) - 0.00048717787299036548*std::cos(19*std::acos(-1 + 2*rpr/Maxx)) + 0.00021535150786820636*std::cos(20*std::acos(-1 + 2*rpr/Maxx)) - 4.7575985736512345e-5*std::cos(21*std::acos(-1 + 2*rpr/Maxx)) - 3.7538066526578624e-5*std::cos(22*std::acos(-1 + 2*rpr/Maxx)) + 6.7304605045776601e-5*std::cos(23*std::acos(-1 + 2*rpr/Maxx)) - 6.5929427728895202e-5*std::cos(24*std::acos(-1 + 2*rpr/Maxx)) + 5.1011723657978163e-5*std::cos(25*std::acos(-1 + 2*rpr/Maxx)) - 3.3360399039120973e-5*std::cos(26*std::acos(-1 + 2*rpr/Maxx)) + 1.8380576200193871e-5*std::cos(27*std::acos(-1 + 2*rpr/Maxx)) - 7.878931024914073e-6*std::cos(28*std::acos(-1 + 2*rpr/Maxx)) + 1.6644334001256316e-6*std::cos(29*std::acos(-1 + 2*rpr/Maxx)) + 1.2960757659209459e-6*std::cos(30*std::acos(-1 + 2*rpr/Maxx)) - 2.178037631388182e-6*std::cos(31*std::acos(-1 + 2*rpr/Maxx)) + 1.9581887788229492e-6*std::cos(32*std::acos(-1 + 2*rpr/Maxx)) - 1.3087487030288237e-6*std::cos(33*std::acos(-1 + 2*rpr/Maxx)) + 6.1900318628474647e-7*std::cos(34*std::acos(-1 + 2*rpr/Maxx)) - 6.8004820184653329e-8*std::cos(35*std::acos(-1 + 2*rpr/Maxx)) - 2.9665542454353737e-7*std::cos(36*std::acos(-1 + 2*rpr/Maxx)) + 4.9510036131430351e-7*std::cos(37*std::acos(-1 + 2*rpr/Maxx)) - 5.7227885640485122e-7*std::cos(38*std::acos(-1 + 2*rpr/Maxx)) + 5.7382153603081516e-7*std::cos(39*std::acos(-1 + 2*rpr/Maxx)) - 5.3555876621860833e-7*std::cos(40*std::acos(-1 + 2*rpr/Maxx)) + 4.8152687055446928e-7*std::cos(41*std::acos(-1 + 2*rpr/Maxx)) - 4.2637502776406691e-7*std::cos(42*std::acos(-1 + 2*rpr/Maxx)) + 3.795025697779657e-7*std::cos(43*std::acos(-1 + 2*rpr/Maxx)) - 3.4898785130852968e-7*std::cos(44*std::acos(-1 + 2*rpr/Maxx)) + 3.4383928358246674e-7*std::cos(45*std::acos(-1 + 2*rpr/Maxx)) - 3.7153750839260132e-7*std::cos(46*std::acos(-1 + 2*rpr/Maxx)) + 4.2910669084185182e-7*std::cos(47*std::acos(-1 + 2*rpr/Maxx)) - 4.9057561744160653e-7*std::cos(48*std::acos(-1 + 2*rpr/Maxx)) - 0.00024358893649518274*std::cos(49*std::acos(-1 + 2*rpr/Maxx)) + 3.14211374001962 - 3.9270934315856834*rpr/Maxx;


         apr = 0.26168769006825793*std::cos(2*std::acos(-1 + 2*rpr/Maxx)) - 0.18455890149887749*std::cos(3*std::acos(-1 + 2*rpr/Maxx)) + 0.10644691856087608*std::cos(4*std::acos(-1 + 2*rpr/Maxx)) - 0.042626922046416882*std::cos(5*std::acos(-1 + 2*rpr/Maxx)) + 0.00056248663002799262*std::cos(6*std::acos(-1 + 2*rpr/Maxx)) + 0.02014476610478377*std::cos(7*std::acos(-1 + 2*rpr/Maxx)) - 0.024957649545127161*std::cos(8*std::acos(-1 + 2*rpr/Maxx)) + 0.020849460083631154*std::cos(9*std::acos(-1 + 2*rpr/Maxx)) - 0.013666027854196879*std::cos(10*std::acos(-1 + 2*rpr/Maxx)) + 0.0070116747178595633*std::cos(11*std::acos(-1 + 2*rpr/Maxx)) - 0.0023727188691233798*std::cos(12*std::acos(-1 + 2*rpr/Maxx)) - 0.00016454816727440303*std::cos(13*std::acos(-1 + 2*rpr/Maxx)) + 0.0011512773593783492*std::cos(14*std::acos(-1 + 2*rpr/Maxx)) - 0.0012464661704612013*std::cos(15*std::acos(-1 + 2*rpr/Maxx)) + 0.00096190009053586376*std::cos(16*std::acos(-1 + 2*rpr/Maxx)) - 0.00060353310259851785*std::cos(17*std::acos(-1 + 2*rpr/Maxx)) + 0.00031014612973111227*std::cos(18*std::acos(-1 + 2*rpr/Maxx)) - 0.00011760298419894008*std::cos(19*std::acos(-1 + 2*rpr/Maxx)) + 1.2972426905246089e-5*std::cos(20*std::acos(-1 + 2*rpr/Maxx)) + 3.1460299297518064e-5*std::cos(21*std::acos(-1 + 2*rpr/Maxx)) - 4.1565148757710908e-5*std::cos(22*std::acos(-1 + 2*rpr/Maxx)) + 3.5796814181932425e-5*std::cos(23*std::acos(-1 + 2*rpr/Maxx)) - 2.5140448974753429e-5*std::cos(24*std::acos(-1 + 2*rpr/Maxx)) + 1.5034014273756846e-5*std::cos(25*std::acos(-1 + 2*rpr/Maxx)) - 7.481343748544104e-6*std::cos(26*std::acos(-1 + 2*rpr/Maxx)) + 2.6858337195837846e-6*std::cos(27*std::acos(-1 + 2*rpr/Maxx)) - 9.9747609096360928e-8*std::cos(28*std::acos(-1 + 2*rpr/Maxx)) - 9.9720849232281744e-7*std::cos(29*std::acos(-1 + 2*rpr/Maxx)) + 1.2325462535254711e-6*std::cos(30*std::acos(-1 + 2*rpr/Maxx)) - 1.0565043162282929e-6*std::cos(31*std::acos(-1 + 2*rpr/Maxx)) + 7.4718662447722397e-7*std::cos(32*std::acos(-1 + 2*rpr/Maxx)) - 4.503170181331782e-7*std::cos(33*std::acos(-1 + 2*rpr/Maxx)) + 2.2444015547850691e-7*std::cos(34*std::acos(-1 + 2*rpr/Maxx)) - 7.8667430662813549e-8*std::cos(35*std::acos(-1 + 2*rpr/Maxx)) - 7.5422423704999825e-10*std::cos(36*std::acos(-1 + 2*rpr/Maxx)) + 3.4257847261451717e-8*std::cos(37*std::acos(-1 + 2*rpr/Maxx)) - 4.0643714033846369e-8*std::cos(38*std::acos(-1 + 2*rpr/Maxx)) + 3.3910631625957225e-8*std::cos(39*std::acos(-1 + 2*rpr/Maxx)) - 2.2969828893536417e-8*std::cos(40*std::acos(-1 + 2*rpr/Maxx)) + 1.26468956304575e-8*std::cos(41*std::acos(-1 + 2*rpr/Maxx)) - 4.9989372907787308e-9*std::cos(42*std::acos(-1 + 2*rpr/Maxx)) + 4.542935982275102e-10*std::cos(43*std::acos(-1 + 2*rpr/Maxx)) + 1.3930406250537502e-9*std::cos(44*std::acos(-1 + 2*rpr/Maxx)) - 1.3113806317392608e-9*std::cos(45*std::acos(-1 + 2*rpr/Maxx)) + 1.7380610306497266e-10*std::cos(46*std::acos(-1 + 2*rpr/Maxx)) + 1.2009518158304569e-9*std::cos(47*std::acos(-1 + 2*rpr/Maxx)) - 2.1730906049510922e-9*std::cos(48*std::acos(-1 + 2*rpr/Maxx)) - 5.8801492099470041e-5*std::cos(49*std::acos(-1 + 2*rpr/Maxx)) + 0.48870556108500007 - 0.63754080718652051*rpr/Maxx;


}
if(rmr < Maxx) 
 { 

         Fmr = -1.0791825914970268*std::cos(2*std::acos(-1 + 2*rmr/Maxx)) + 0.80285889978116265*std::cos(3*std::acos(-1 + 2*rmr/Maxx)) - 0.54006781332116993*std::cos(4*std::acos(-1 + 2*rmr/Maxx)) + 0.33041822203594301*std::cos(5*std::acos(-1 + 2*rmr/Maxx)) - 0.18265025176158761*std::cos(6*std::acos(-1 + 2*rmr/Maxx)) + 0.088143335337428227*std::cos(7*std::acos(-1 + 2*rmr/Maxx)) - 0.032914921489218332*std::cos(8*std::acos(-1 + 2*rmr/Maxx)) + 0.0039265669028494682*std::cos(9*std::acos(-1 + 2*rmr/Maxx)) + 0.0088684564531565507*std::cos(10*std::acos(-1 + 2*rmr/Maxx)) - 0.012487917176595156*std::cos(11*std::acos(-1 + 2*rmr/Maxx)) + 0.011512829661162502*std::cos(12*std::acos(-1 + 2*rmr/Maxx)) - 0.008736996468339267*std::cos(13*std::acos(-1 + 2*rmr/Maxx)) + 0.0057139888532271449*std::cos(14*std::acos(-1 + 2*rmr/Maxx)) - 0.0031811020490695868*std::cos(15*std::acos(-1 + 2*rmr/Maxx)) + 0.0013766142424895406*std::cos(16*std::acos(-1 + 2*rmr/Maxx)) - 0.0002666432511588962*std::cos(17*std::acos(-1 + 2*rmr/Maxx)) - 0.00030193860802503469*std::cos(18*std::acos(-1 + 2*rmr/Maxx)) + 0.00050814209550134604*std::cos(19*std::acos(-1 + 2*rmr/Maxx)) - 0.00050768215750462975*std::cos(20*std::acos(-1 + 2*rmr/Maxx)) + 0.00041424714553496556*std::cos(21*std::acos(-1 + 2*rmr/Maxx)) - 0.00029887944196967109*std::cos(22*std::acos(-1 + 2*rmr/Maxx)) + 0.00019832092411475877*std::cos(23*std::acos(-1 + 2*rmr/Maxx)) - 0.00012591855763192938*std::cos(24*std::acos(-1 + 2*rmr/Maxx)) + 8.1570792144811111e-5*std::cos(25*std::acos(-1 + 2*rmr/Maxx)) - 5.9127074289927473e-5*std::cos(26*std::acos(-1 + 2*rmr/Maxx)) + 5.1061133080856913e-5*std::cos(27*std::acos(-1 + 2*rmr/Maxx)) - 5.0876555611187313e-5*std::cos(28*std::acos(-1 + 2*rmr/Maxx)) + 5.396483258066076e-5*std::cos(29*std::acos(-1 + 2*rmr/Maxx)) - 5.7571240465195683e-5*std::cos(30*std::acos(-1 + 2*rmr/Maxx)) + 6.0370580947335939e-5*std::cos(31*std::acos(-1 + 2*rmr/Maxx)) - 6.1961389745610052e-5*std::cos(32*std::acos(-1 + 2*rmr/Maxx)) + 6.2437584665536727e-5*std::cos(33*std::acos(-1 + 2*rmr/Maxx)) - 6.2089489739611259e-5*std::cos(34*std::acos(-1 + 2*rmr/Maxx)) + 6.1228909150320976e-5*std::cos(35*std::acos(-1 + 2*rmr/Maxx)) - 6.0107686167987229e-5*std::cos(36*std::acos(-1 + 2*rmr/Maxx)) + 5.8895261093941387e-5*std::cos(37*std::acos(-1 + 2*rmr/Maxx)) - 5.7687291804345286e-5*std::cos(38*std::acos(-1 + 2*rmr/Maxx)) + 5.6525168130025505e-5*std::cos(39*std::acos(-1 + 2*rmr/Maxx)) - 5.5417299469798231e-5*std::cos(40*std::acos(-1 + 2*rmr/Maxx)) + 5.4353656013819696e-5*std::cos(41*std::acos(-1 + 2*rmr/Maxx)) - 5.3320161792283494e-5*std::cos(42*std::acos(-1 + 2*rmr/Maxx)) + 5.2298314924745237e-5*std::cos(43*std::acos(-1 + 2*rmr/Maxx)) - 5.1285003905179871e-5*std::cos(44*std::acos(-1 + 2*rmr/Maxx)) + 5.0263916858819813e-5*std::cos(45*std::acos(-1 + 2*rmr/Maxx)) - 4.9283781304740315e-5*std::cos(46*std::acos(-1 + 2*rmr/Maxx)) + 4.8307125567192966e-5*std::cos(47*std::acos(-1 + 2*rmr/Maxx)) - 4.7592632780077415e-5*std::cos(48*std::acos(-1 + 2*rmr/Maxx)) + 2.3453717968781318e-5*std::cos(49*std::acos(-1 + 2*rmr/Maxx)) + 3.0117542335016725 + 2.5949803874590023*rmr/Maxx;


         Chimr = 1.1795075506129167*std::cos(2*std::acos(-1 + 2*rmr/Maxx)) - 0.53917011624401534*std::cos(3*std::acos(-1 + 2*rmr/Maxx)) + 0.16962696949652423*std::cos(4*std::acos(-1 + 2*rmr/Maxx)) - 0.0037342126895790538*std::cos(5*std::acos(-1 + 2*rmr/Maxx)) - 0.053759834803502007*std::cos(6*std::acos(-1 + 2*rmr/Maxx)) + 0.062073415559850888*std::cos(7*std::acos(-1 + 2*rmr/Maxx)) - 0.051729276271027136*std::cos(8*std::acos(-1 + 2*rmr/Maxx)) + 0.036745573212263516*std::cos(9*std::acos(-1 + 2*rmr/Maxx)) - 0.023006732885378415*std::cos(10*std::acos(-1 + 2*rmr/Maxx)) + 0.012502022934543887*std::cos(11*std::acos(-1 + 2*rmr/Maxx)) - 0.0054260815940578196*std::cos(12*std::acos(-1 + 2*rmr/Maxx)) + 0.0012235286986818584*std::cos(13*std::acos(-1 + 2*rmr/Maxx)) + 0.00088379036280092978*std::cos(14*std::acos(-1 + 2*rmr/Maxx)) - 0.0016379461520883594*std::cos(15*std::acos(-1 + 2*rmr/Maxx)) + 0.0016334052552107966*std::cos(16*std::acos(-1 + 2*rmr/Maxx)) - 0.0012885284225685094*std::cos(17*std::acos(-1 + 2*rmr/Maxx)) + 0.00086167705025049965*std::cos(18*std::acos(-1 + 2*rmr/Maxx)) - 0.00048717787299036548*std::cos(19*std::acos(-1 + 2*rmr/Maxx)) + 0.00021535150786820636*std::cos(20*std::acos(-1 + 2*rmr/Maxx)) - 4.7575985736512345e-5*std::cos(21*std::acos(-1 + 2*rmr/Maxx)) - 3.7538066526578624e-5*std::cos(22*std::acos(-1 + 2*rmr/Maxx)) + 6.7304605045776601e-5*std::cos(23*std::acos(-1 + 2*rmr/Maxx)) - 6.5929427728895202e-5*std::cos(24*std::acos(-1 + 2*rmr/Maxx)) + 5.1011723657978163e-5*std::cos(25*std::acos(-1 + 2*rmr/Maxx)) - 3.3360399039120973e-5*std::cos(26*std::acos(-1 + 2*rmr/Maxx)) + 1.8380576200193871e-5*std::cos(27*std::acos(-1 + 2*rmr/Maxx)) - 7.878931024914073e-6*std::cos(28*std::acos(-1 + 2*rmr/Maxx)) + 1.6644334001256316e-6*std::cos(29*std::acos(-1 + 2*rmr/Maxx)) + 1.2960757659209459e-6*std::cos(30*std::acos(-1 + 2*rmr/Maxx)) - 2.178037631388182e-6*std::cos(31*std::acos(-1 + 2*rmr/Maxx)) + 1.9581887788229492e-6*std::cos(32*std::acos(-1 + 2*rmr/Maxx)) - 1.3087487030288237e-6*std::cos(33*std::acos(-1 + 2*rmr/Maxx)) + 6.1900318628474647e-7*std::cos(34*std::acos(-1 + 2*rmr/Maxx)) - 6.8004820184653329e-8*std::cos(35*std::acos(-1 + 2*rmr/Maxx)) - 2.9665542454353737e-7*std::cos(36*std::acos(-1 + 2*rmr/Maxx)) + 4.9510036131430351e-7*std::cos(37*std::acos(-1 + 2*rmr/Maxx)) - 5.7227885640485122e-7*std::cos(38*std::acos(-1 + 2*rmr/Maxx)) + 5.7382153603081516e-7*std::cos(39*std::acos(-1 + 2*rmr/Maxx)) - 5.3555876621860833e-7*std::cos(40*std::acos(-1 + 2*rmr/Maxx)) + 4.8152687055446928e-7*std::cos(41*std::acos(-1 + 2*rmr/Maxx)) - 4.2637502776406691e-7*std::cos(42*std::acos(-1 + 2*rmr/Maxx)) + 3.795025697779657e-7*std::cos(43*std::acos(-1 + 2*rmr/Maxx)) - 3.4898785130852968e-7*std::cos(44*std::acos(-1 + 2*rmr/Maxx)) + 3.4383928358246674e-7*std::cos(45*std::acos(-1 + 2*rmr/Maxx)) - 3.7153750839260132e-7*std::cos(46*std::acos(-1 + 2*rmr/Maxx)) + 4.2910669084185182e-7*std::cos(47*std::acos(-1 + 2*rmr/Maxx)) - 4.9057561744160653e-7*std::cos(48*std::acos(-1 + 2*rmr/Maxx)) - 0.00024358893649518274*std::cos(49*std::acos(-1 + 2*rmr/Maxx)) + 3.14211374001962 - 3.9270934315856834*rmr/Maxx;


         amr = 0.26168769006825793*std::cos(2*std::acos(-1 + 2*rmr/Maxx)) - 0.18455890149887749*std::cos(3*std::acos(-1 + 2*rmr/Maxx)) + 0.10644691856087608*std::cos(4*std::acos(-1 + 2*rmr/Maxx)) - 0.042626922046416882*std::cos(5*std::acos(-1 + 2*rmr/Maxx)) + 0.00056248663002799262*std::cos(6*std::acos(-1 + 2*rmr/Maxx)) + 0.02014476610478377*std::cos(7*std::acos(-1 + 2*rmr/Maxx)) - 0.024957649545127161*std::cos(8*std::acos(-1 + 2*rmr/Maxx)) + 0.020849460083631154*std::cos(9*std::acos(-1 + 2*rmr/Maxx)) - 0.013666027854196879*std::cos(10*std::acos(-1 + 2*rmr/Maxx)) + 0.0070116747178595633*std::cos(11*std::acos(-1 + 2*rmr/Maxx)) - 0.0023727188691233798*std::cos(12*std::acos(-1 + 2*rmr/Maxx)) - 0.00016454816727440303*std::cos(13*std::acos(-1 + 2*rmr/Maxx)) + 0.0011512773593783492*std::cos(14*std::acos(-1 + 2*rmr/Maxx)) - 0.0012464661704612013*std::cos(15*std::acos(-1 + 2*rmr/Maxx)) + 0.00096190009053586376*std::cos(16*std::acos(-1 + 2*rmr/Maxx)) - 0.00060353310259851785*std::cos(17*std::acos(-1 + 2*rmr/Maxx)) + 0.00031014612973111227*std::cos(18*std::acos(-1 + 2*rmr/Maxx)) - 0.00011760298419894008*std::cos(19*std::acos(-1 + 2*rmr/Maxx)) + 1.2972426905246089e-5*std::cos(20*std::acos(-1 + 2*rmr/Maxx)) + 3.1460299297518064e-5*std::cos(21*std::acos(-1 + 2*rmr/Maxx)) - 4.1565148757710908e-5*std::cos(22*std::acos(-1 + 2*rmr/Maxx)) + 3.5796814181932425e-5*std::cos(23*std::acos(-1 + 2*rmr/Maxx)) - 2.5140448974753429e-5*std::cos(24*std::acos(-1 + 2*rmr/Maxx)) + 1.5034014273756846e-5*std::cos(25*std::acos(-1 + 2*rmr/Maxx)) - 7.481343748544104e-6*std::cos(26*std::acos(-1 + 2*rmr/Maxx)) + 2.6858337195837846e-6*std::cos(27*std::acos(-1 + 2*rmr/Maxx)) - 9.9747609096360928e-8*std::cos(28*std::acos(-1 + 2*rmr/Maxx)) - 9.9720849232281744e-7*std::cos(29*std::acos(-1 + 2*rmr/Maxx)) + 1.2325462535254711e-6*std::cos(30*std::acos(-1 + 2*rmr/Maxx)) - 1.0565043162282929e-6*std::cos(31*std::acos(-1 + 2*rmr/Maxx)) + 7.4718662447722397e-7*std::cos(32*std::acos(-1 + 2*rmr/Maxx)) - 4.503170181331782e-7*std::cos(33*std::acos(-1 + 2*rmr/Maxx)) + 2.2444015547850691e-7*std::cos(34*std::acos(-1 + 2*rmr/Maxx)) - 7.8667430662813549e-8*std::cos(35*std::acos(-1 + 2*rmr/Maxx)) - 7.5422423704999825e-10*std::cos(36*std::acos(-1 + 2*rmr/Maxx)) + 3.4257847261451717e-8*std::cos(37*std::acos(-1 + 2*rmr/Maxx)) - 4.0643714033846369e-8*std::cos(38*std::acos(-1 + 2*rmr/Maxx)) + 3.3910631625957225e-8*std::cos(39*std::acos(-1 + 2*rmr/Maxx)) - 2.2969828893536417e-8*std::cos(40*std::acos(-1 + 2*rmr/Maxx)) + 1.26468956304575e-8*std::cos(41*std::acos(-1 + 2*rmr/Maxx)) - 4.9989372907787308e-9*std::cos(42*std::acos(-1 + 2*rmr/Maxx)) + 4.542935982275102e-10*std::cos(43*std::acos(-1 + 2*rmr/Maxx)) + 1.3930406250537502e-9*std::cos(44*std::acos(-1 + 2*rmr/Maxx)) - 1.3113806317392608e-9*std::cos(45*std::acos(-1 + 2*rmr/Maxx)) + 1.7380610306497266e-10*std::cos(46*std::acos(-1 + 2*rmr/Maxx)) + 1.2009518158304569e-9*std::cos(47*std::acos(-1 + 2*rmr/Maxx)) - 2.1730906049510922e-9*std::cos(48*std::acos(-1 + 2*rmr/Maxx)) - 5.8801492099470041e-5*std::cos(49*std::acos(-1 + 2*rmr/Maxx)) + 0.48870556108500007 - 0.63754080718652051*rmr/Maxx;


}
        
        amrex::Real Fplmod = Fpl;

        amrex::Real Fmlmod = Fml;

        amrex::Real Fprmod = Fpr;

        amrex::Real Fmrmod = Fmr;

        amrex::Real Chiplmod = Chipl;

        amrex::Real Chimlmod = Chiml;

        amrex::Real Chiprmod = Chipr;

        amrex::Real Chimrmod = Chimr;
        
        
        
        
        amrex::Real A_L_0pl = 1/e*(ylprime/(rpl*rpl+0.0001))*(1-apl);
        
        amrex::Real A_L_1pl = -1/e*(xpl/(rpl*rpl+0.0001))*(1-apl);
        
        amrex::Real A_L_2pl = 0;
        
        amrex::Real A_L_0ml = 1/e*(ylprime/(rml*rml+0.0001))*(1-aml);
        
        amrex::Real A_L_1ml = -1/e*(xml/(rml*rml+0.0001))*(1-aml);
        
        amrex::Real A_L_2ml = 0;
        
        amrex::Real A_L_0pr = 1/e*(yrprime/(rpr*rpr+0.0001))*(1-apr);
        
        amrex::Real A_L_1pr = -1/e*(xpr/(rpr*rpr+0.0001))*(1-apr);
        
        amrex::Real A_L_2pr = 0;
        
        amrex::Real A_L_0mr = 1/e*(yrprime/(rmr*rmr+0.0001))*(1-amr);
        
        amrex::Real A_L_1mr = -1/e*(xmr/(rmr*rmr+0.0001))*(1-amr);
        
        amrex::Real A_L_2mr = 0;
        
        /*
        amrex::Real Fplmod = Fpl*(rpl < Maxx) + (rpl >= Maxx);

        amrex::Real Fmlmod = Fml*(rml < Maxx) + (rml >= Maxx);

        amrex::Real Fprmod = Fpr*(rpr < Maxx) + (rpr >= Maxx);

        amrex::Real Fmrmod = Fmr*(rmr < Maxx) + (rmr >= Maxx);

        amrex::Real Chiplmod = Chipl*(rpl < Maxx);

        amrex::Real Chimlmod = Chiml*(rml < Maxx);

        amrex::Real Chiprmod = Chipr*(rpr < Maxx);

        amrex::Real Chimrmod = Chimr*(rmr < Maxx);
        */
        /*
        amrex::Real Phirealp = Fplmod*Fprmod*(xpl*xpr - ylprime*yrprime)/(rpl*rpr*v);

        amrex::Real Phiimaginaryp = Fplmod*Fprmod*(xpl*yrprime + xpr*ylprime)/(rpl*rpr*v);

        amrex::Real Phirealm = Fmlmod*Fmrmod*(xml*xmr - ylprime*yrprime)/(rml*rmr*v);

        amrex::Real Phiimaginarym = Fmlmod*Fmrmod*(xml*yrprime + xmr*ylprime)/(rml*rmr*v);

        amrex::Real Chip = Chiplmod + Chiprmod;

        amrex::Real Chim = Chimlmod + Chimrmod;
        
        amrex::Real ATp = gammaboost*vboost*A_L_0pl - gammaboost*vboost*A_L_0pr;
        
        amrex::Real A_L_0p  = gammaboost*A_L_0pl + gammaboost*A_L_0pr;
        
        amrex::Real A_L_1p  = A_L_1pl + A_L_1pr;
        
        amrex::Real A_L_2p  = A_L_2pl + A_L_2pr;
        
        amrex::Real ATm = gammaboost*vboost*A_L_0ml - gammaboost*vboost*A_L_0mr;
        
        amrex::Real A_L_0m  = gammaboost*A_L_0ml + gammaboost*A_L_0mr;
        
        amrex::Real A_L_1m  = A_L_1ml + A_L_1mr;
        
        amrex::Real A_L_2m  = A_L_2ml + A_L_2mr;
        */
        
        amrex::Real Phirealp = Fplmod*Fprmod*(xpl*xpr - ylprime*yrprime)/(rpl*rpr*v);

        amrex::Real Phiimaginaryp = Fplmod*Fprmod*(xpl*yrprime + xpr*ylprime)/(rpl*rpr*v);

        amrex::Real Phirealm = Fmlmod*Fmrmod*(xml*xmr - ylprime*yrprime)/(rml*rmr*v);

        amrex::Real Phiimaginarym = Fmlmod*Fmrmod*(xml*yrprime + xmr*ylprime)/(rml*rmr*v);

        amrex::Real Chip = Chiplmod - Chiprmod;
        
        amrex::Real Chim = Chimlmod - Chimrmod;
        
        amrex::Real ATp = -gammaboost*vboost*A_L_0pl + gammaboost*vboost*A_L_0pr;
        
        amrex::Real A_L_0p  = gammaboost*A_L_0pl + gammaboost*A_L_0pr;
        
        amrex::Real A_L_1p  = std::cos(angle1)*A_L_1pl + std::sin(angle1)*A_L_2pl + std::cos(angle2)*A_L_1pr + std::sin(angle2)*A_L_2pr;
        
        amrex::Real A_L_2p  = std::cos(angle1)*A_L_2pl - std::sin(angle1)*A_L_1pl + std::cos(angle2)*A_L_2pr - std::sin(angle2)*A_L_1pr;
        
        amrex::Real ATm = -gammaboost*vboost*A_L_0ml + gammaboost*vboost*A_L_0mr;
        
        amrex::Real A_L_0m  = gammaboost*A_L_0ml + gammaboost*A_L_0mr;
       
        amrex::Real A_L_1m  = std::cos(angle1)*A_L_1ml + std::sin(angle1)*A_L_2ml + std::cos(angle2)*A_L_1mr + std::sin(angle2)*A_L_2mr;
        
        amrex::Real A_L_2m  = std::cos(angle1)*A_L_2ml - std::sin(angle1)*A_L_1ml + std::cos(angle2)*A_L_2mr - std::sin(angle2)*A_L_1mr;
        
        
        /*
        amrex::Real Phirealp = Fplmod*Fprmod*(xpl*xpr - y*y)/(rpl*rpr*v);

        amrex::Real Phiimaginaryp = Fplmod*Fprmod*(xpl*y + xpr*y)/(rpl*rpr*v);

        amrex::Real Phirealm = Fmlmod*Fmrmod*(xml*xmr - y*y)/(rml*rmr*v);

        amrex::Real Phiimaginarym = Fmlmod*Fmrmod*(xml*y + xmr*y)/(rml*rmr*v);

        amrex::Real Chip = Chiplmod + Chiprmod;

        amrex::Real Chim = Chimlmod + Chimrmod;
        */
        
        state_fab(i, j, k, Idx::PhiReal) = 0.5*Phirealm + 0.5*Phirealp;

        state_fab(i, j, k, Idx::PhiImaginary) = 0.5*Phiimaginarym + 0.5*Phiimaginaryp;

        state_fab(i, j, k, Idx::Chi) = 0.5*Chim + 0.5*Chip;
        
        state_fab(i, j, k, Idx::Chi0) = 0;
        
        state_fab(i, j, k, Idx::Chi1) = 0;
        
        state_fab(i, j, k, Idx::AT) = 0.5*ATm + 0.5*ATp;
        
        state_fab(i, j, k, Idx::A_L_0) = 0.5*A_L_0m + 0.5*A_L_0p;
        
        state_fab(i, j, k, Idx::A_L_1) = 0.5*A_L_1m + 0.5*A_L_1p;
        
        state_fab(i, j, k, Idx::A_L_2) = 0.5*A_L_2m + 0.5*A_L_2p;
        
        state_fab(i, j, k, Idx::PiPhiReal) = (-Phirealm + Phirealp)/(-tm + tp);

        state_fab(i, j, k, Idx::PiPhiImaginary) = (-Phiimaginarym + Phiimaginaryp)/(-tm + tp);

        state_fab(i, j, k, Idx::PiChi) = (-Chim + Chip)/(-tm + tp);
        
        state_fab(i, j, k, Idx::PiChi0) = 0;
        
        state_fab(i, j, k, Idx::PiChi1) = 0;
        
        state_fab(i, j, k, Idx::PiAT) = (-ATm + ATp)/(-tm + tp);
        
        state_fab(i, j, k, Idx::PiA_L_0) = (-A_L_0m + A_L_0p)/(-tm + tp);
        
        state_fab(i, j, k, Idx::PiA_L_1) = (-A_L_1m + A_L_1p)/(-tm + tp);
        
        state_fab(i, j, k, Idx::PiA_L_2) = (-A_L_2m + A_L_2p)/(-tm + tp);
        
}
#endif