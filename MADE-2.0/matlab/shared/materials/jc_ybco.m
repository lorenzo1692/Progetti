function [jc, Birr] = jc_ybco(B, T, opt)
%JC_YBCO REBCO tape critical current density fit (A/m^2).
%
%   [jc, Birr] = JC_YBCO(B, T, opt) with opt a named dataset ('so', 'sp',
%   'sp2', 'sst', 'sst2', 'sst_parr', 'kit') or an explicit
%   [A alpha beta p q Birr0 Tc sf] vector.
%
%   Relocated unchanged from the CS legacy archive (jc_ybco.m). Called by
%   Ic_sst33.m.

if nargin < 3
    opt = 'so';
end

if ischar(opt)
    switch opt
        case 'so' % ENEA
            A=1.8526e12; alpha=1.4264; betta=2.3795; p=0.39042; q=1.0426; Birr0=119.96; Tc=92.091; sf=7.8992e-3; %#ok<NASGU>
        case 'sp' % Fermilab'09
            A=2.0441e12; alpha=1.5254; betta=2.1040; p=0.4616; q=0.4134; Birr0=120.02; Tc=92.90; sf=10.144e-3; %#ok<NASGU>
        case 'sp2' % Rainer, low temperature
            A=2.0381e12; alpha=1.5922; betta=2.0621; p=0.47257; q=1.3174; Birr0=119.97; Tc=92.914; sf=36.737e-3; %#ok<NASGU>
        case 'sst' % from brochure
            A=3.5739e12; alpha=1.5181; betta=2.3251; p=0.50215; q=1.6983; Birr0=120; Tc=92.833; sf=3.4675e-3; %#ok<NASGU>
        case 'sst2' % from brochure
            A=3.5739e12*1.0611; alpha=1.5181; betta=2.3251; p=0.50215; q=1.6983; Birr0=120; Tc=92.833; sf=3.4675e-3; %#ok<NASGU>
        case 'sst_parr' % 4.2K, 8-15T fit
            A=1.6096e13; alpha=0.9181; betta=3.4251; p=0.60215; q=0.7983; Birr0=120; Tc=92.833; sf=3.4675e-3; %#ok<NASGU>
        case 'kit' % Reinhard option
            A=2.7310e12; alpha=1.54121; betta=1.96679; p=0.5875; q=1.7; Birr0=132.5; Tc=90; sf=5e-3; %#ok<NASGU>
        otherwise
            error('jc_ybco:unknown_dataset', 'Unknown dataset option: %s', opt);
    end
else
    A=opt(1); alpha=opt(2); betta=opt(3); p=opt(4); q=opt(5); Birr0=opt(6); Tc=opt(7); sf=opt(8); %#ok<NASGU>
end

t = T/Tc;
Birr = Birr0*(1-t).^alpha;
Birr(t>1) = 1e-16;
B(B<sf) = sf;
b = B./Birr;

jc = A./B.*((Birr/Birr0).^betta).*(b.^p).*(1-b).^q;
jc(b>1) = 1e-16;
end
