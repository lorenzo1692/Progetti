function d = tf3d_design_geometry(row,p)
%TF3D_DESIGN_GEOMETRY Real turns, in metres. Coil 1 lies in the YZ plane.
% X_3D = r*e_r + v*e_phi + z*e_z; at phi=pi/2, x_2D = v = -X_3D.
% No machine-radius substitution is made for the selected row.Ri_.
if ~isstruct(row), row=table2struct(row); end
if numel(row)~=1, error('tf3d:row','Select exactly one design.'); end
required={'Iop','n_layers','n_turns','Cond_w','Cond_h','JT','Ri_','Rk_'};
for i=1:numel(required)
    if ~isfield(row,required{i}), error('tf3d:row','Missing row.%s',required{i}); end
end
nl=row.n_layers;
if ~isscalar(nl)||~isfinite(nl)||nl<1||nl~=round(nl), error('tf3d:geometry','Invalid n_layers.'); end
for key={'n_turns','Cond_w','Cond_h','JT'}
    v=row.(key{1});
    if ~isnumeric(v)||numel(v)<nl||any(~isfinite(v(1:nl)))||any(v(1:nl)<=0)
        error('tf3d:geometry','Invalid %s.',key{1});
    end
end
nt=row.n_turns(1:nl);
if any(nt~=round(nt)), error('tf3d:geometry','Turn counts must be integers.'); end
if ~isscalar(row.Iop)||~isfinite(row.Iop)||row.Iop<=0, error('tf3d:geometry','Iop must be positive.'); end
if p.n_TF<2||p.n_TF~=round(p.n_TF), error('tf3d:geometry','n_TF must be an integer >=2.'); end
if p.GoundIns<0||p.INS_grades<0||p.dr_plasma_side<0||p.Increm<=0||p.turn_insulation_nominal<0
    error('tf3d:geometry','Invalid insulation or case dimensions.');
end
tins=p.turn_insulation_nominal*p.Increm;
% Conductor shape from the SOLUTION (RIS stays RIS, see WP_TURN_GEOMETRY).
tg=wp_turn_geometry(row,p);
d.shape_name=tg.shape_name;d.cable_is_round=tg.is_round;
% Gauss sources: a square of side s has second moment s^2/12 per axis, a
% circle of diameter D has D^2/16, so a RIS cable is represented by the
% equivalent square of side D*sqrt(3)/2 (same centroidal second moment).
qf=1-(1-sqrt(3)/2)*tg.is_round;
d.qw=[];d.qh=[];
d.x=[];d.y=[];d.layer=[];d.w=[];d.h=[];d.cell_h=[];d.cell_w=[];
Re=row.Ri_-p.dr_plasma_side-p.GoundIns;
d.wp_plasma_radius=Re;
for k=1:nl
    cw=row.Cond_w(k);ch=row.Cond_h(k);j=row.JT(k);
    d.x=[d.x ((1:nt(k))-(nt(k)+1)/2)*cw]; %#ok<AGROW>
    d.y=[d.y repmat(Re-ch/2,1,nt(k))]; %#ok<AGROW>
    d.layer=[d.layer repmat(k,1,nt(k))]; %#ok<AGROW>
    d.w=[d.w repmat(cw-2*j-2*tins,1,nt(k))]; %#ok<AGROW>
    d.h=[d.h repmat(ch-2*j-2*tins,1,nt(k))]; %#ok<AGROW>
    d.qw=[d.qw repmat(qf*tg.cab_w(k),1,nt(k))]; %#ok<AGROW>
    d.qh=[d.qh repmat(qf*tg.cab_h(k),1,nt(k))]; %#ok<AGROW>
    d.cell_h=[d.cell_h repmat(ch,1,nt(k))]; %#ok<AGROW>
    d.cell_w=[d.cell_w repmat(cw,1,nt(k))]; %#ok<AGROW>
    if k<nl, Re=Re-ch-p.INS_grades; else, Re=Re-ch; end
end
if any(d.w<=0|d.h<=0), error('tf3d:geometry','Non-positive cable cross-section.'); end
d.ground_insulation=p.GoundIns;d.wp_back_radius=Re; d.r_c=mean(d.y);d.d_c=row.Ri_-d.r_c;
g=compute_operating_params(p);
d.r1=d.r_c;d.r2=g.RTFo+d.d_c;
d.n_turns=sum(nt);d.Iop=row.Iop;d.NI=d.Iop*d.n_turns;
d.u=d.r_c-d.y;d.v=d.x;
d.wp_u_min=d.r_c-d.wp_plasma_radius;
d.wp_u_max=d.r_c-d.wp_back_radius;
d.case_u_min=-d.d_c;
% Nose arc radius at the centre plane, not the wedge-face radius Rk_.
d.case_inner_radius=row.Rk_/cos(pi/p.n_TF);
d.case_u_max=d.r_c-d.case_inner_radius;
d.nose_il=d.wp_back_radius-p.GoundIns-d.case_inner_radius;
d.wp_width=max(nt(:)'.*reshape(row.Cond_w(1:nl),1,[]));
d.wp_depth=d.wp_plasma_radius-d.wp_back_radius;
d.machine_inner_radius=g.RTFi;d.row_inner_radius=row.Ri_;
if any(~isfinite([d.r1 d.r2 d.nose_il]))||d.r1<=0||d.r2<=d.r1||d.nose_il<0
    error('tf3d:geometry','Invalid current-line radii or negative nose thickness.');
end
end
