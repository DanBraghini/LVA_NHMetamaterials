function  function_unfolded_dispersion_diagrams(pb,kLv,w_real,w_imag,varargin)
% function_unfolded_dispersion_diagrams(pb,kLv,w_real_TM,w_imag_TM,w_real_TM_passive,w_imag_TM_passive,varargin)
% this function plots dispersion diagram of a non-Hermitian system for up to
% ten bulk bands
% last update: 08, Janyary, 2022
%%
if nargin > 4
    options = struct(varargin{:});
else
    options = [];
end

if ~isfield(options,'setplot')
    options.setplot = [1 1 1 0];
end
if ~isfield(options,'periodic')
    options.periodic = 0;
end
if ~isfield(options,'plotpassive')
    options.plotpassive = 0;
else
    if isfield(options,'w_real_passive') && isfield(options,'w_imag_passive')
        w_real_passive = options.w_real_passive;
        w_imag_passive = options.w_imag_passive;
    end
end
if ~isfield(options,'plotcolors')
    options.plotcolors = 0;
end

%%
% Indexing the boundaries between the Brillouin Zones
dk = kLv(2)-kLv(1);
kind=zeros(1,2*pb);
j=1;
for i=1:2:2*pb
   kindp=find(abs(kLv + j*pi) <= dk/2);
   if length(kindp) ~=1
       kind(i)=kindp(1);
   else
       kind(i)=kindp;
   end
   kindm=find(abs(kLv - j*pi) <= dk/2);
   if length(kindm) ~=1
       kind(i+1)=kindm(1);
   else
       kind(i+1)=kindm;
   end
   j=j+1;
end

kv=kLv/pi;
if options.setplot(1)==1
    % plot frequency by normalized wavenumber 
    %each PB is uniquely associated with one BZ
    figure
    plot(kv(kind(1):kind(2)),w_real{1}(kind(1):kind(2)),'k')
    hold on
    if options.plotpassive == 1
        plot(kv(kind(1):kind(2)),w_real_passive{1}(kind(1):kind(2)),'r')
    end
    if options.periodic~=0    
        for i=1:pb
            plot(kv, w_real{i},'--k');
            if options.plotpassive == 1
                plot(kv, w_real_passive{i},'--r');
            end
        end
    end
    for i=2:pb
        plot(kv(kind(2*i-1):kind(2*i-3)),w_real{i}(kind(2*i-1):kind(2*i-3)),'k',...
             kv(kind(2*i-2):kind(2*i)),w_real{i}(kind(2*i-2):kind(2*i)),'k')
        if options.plotpassive == 1
        plot(kv(kind(2*i-1):kind(2*i-3)),w_real_passive{i}(kind(2*i-1):kind(2*i-3)),'r',...
             kv(kind(2*i-2):kind(2*i)),w_real_passive{i}(kind(2*i-2):kind(2*i)),'r')
        end
    end
    ylabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$k L_c/ \pi$', 'interpreter', 'latex', 'fontsize', 15)
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    box on
    set(gcf, 'Color', 'w');
    ylim([0 w_real{pb}(end)])
    xlim([-pb pb])
    grid on
    % export_fig fig2_unfolded1.pdf
end
if options.setplot(2)==1
    figure
    % plot damping_factor by normalized wavenumber 
    %each PB is uniquely associated with one BZ
    plot(kv(kind(1):kind(2)),w_imag{1}(kind(1):kind(2)),'k')
     hold on
     if options.plotpassive == 1
        plot(kv(kind(1):kind(2)),w_imag_passive{1}(kind(1):kind(2)),'r')
     end
    for i=2:pb
        plot(kv(kind(2*i-1):kind(2*i-3)),w_imag{i}(kind(2*i-1):kind(2*i-3)),'k',...
             kv(kind(2*i-2):kind(2*i)),w_imag{i}(kind(2*i-2):kind(2*i)),'k')
     if options.plotpassive == 1
          plot(kv(kind(2*i-1):kind(2*i-3)),w_imag_passive{i}(kind(2*i-1):kind(2*i-3)),'r',...
             kv(kind(2*i-2):kind(2*i)),w_imag_passive{i}(kind(2*i-2):kind(2*i)),'r')
     end
    end

    xlim([-pb pb])
    ylabel('$\Im(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$k L_c/ \pi$', 'interpreter', 'latex', 'fontsize', 15)
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    box on
    grid on
    set(gcf, 'Color', 'w');
    % export_fig fig2_unfolded2.pdf
end
%% Plotting complex plane
if options.setplot(3)==1
    figure
    for j=1:pb
        if options.plotcolors == 1
            scatter(w_imag{j}, w_real{j})
        else
            plot(w_imag{j}, w_real{j},'k')
        end
        hold on
        if options.plotpassive == 1
            plot(w_imag_passive{j}, w_real_passive{j},'r')
        end
     end
    ylabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$\Im(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    box on
    grid on
    set(gcf, 'Color', 'w');
    % export_fig fig2_unfolded3.pdf
end
%% Plotting 3D reciprocal space
if options.setplot(4)==1
    figure
     plot3(w_real{1}, w_imag{1},kv,'r','MarkerSize',8,'LineWidth',2)
     hold on
      plot3(w_real{2}, w_imag{2},kv,'b','MarkerSize',8,'LineWidth',2)
      plot3(w_real{3}, w_imag{3},kv,'c','MarkerSize',8,'LineWidth',2)
      plot3(w_real{4}, w_imag{4},kv,'m','MarkerSize',8,'LineWidth',2)
    for j=1:pb
    %    plot3(w_real{j}, w_imag{j},kv,'b','MarkerSize',8,'LineWidth',2)
   %     hold on
%        scatter3(w_real{j}, w_imag{j},kv,'m')
        plot(w_real{j}, w_imag{j},'k','MarkerSize',8,'LineWidth',2)
        scatter(w_real{j}, w_imag{j},'g')
    end
    xlabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$\Im(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 20)
    zlabel('$k L_c/ \pi$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca,'TickLabelInterpreter','Latex','fontsize',20);
    box on
    set(gcf, 'Color', 'w');
    grid on
    
end

end