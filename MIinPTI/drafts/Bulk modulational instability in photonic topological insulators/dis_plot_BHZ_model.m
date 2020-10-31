clear all
close all;
%%%%%%%%BULK STATE (NonDis)

%���������
% k_min=-1;
% delta_k=0.15;
% delta_w=delta_k;
% k_max=1;
k_min=-2;
delta_k=0.04;
delta_w=delta_k;
k_max=2;
beta=  1;
% I_0=1;
I_0=3;
M=-1;
g=1;


%����������� ������� w_solved(������ �� �,������������� ����� ������ ������� ��������� (��� �� ���������))
k=[k_min:delta_k:k_max];
syms w
for i=1:1:length(k)
    eq=(w+g*I_0/2).^2 - k(i).^2 - ((M+beta*k(i).^2)-g*I_0*(M+beta*k(i).^2)./2./(w+g*I_0)).^2==0;
    w_solved_str=solve(eq,w);
    %%%%%
    %%%%%%
    %%%%%%% ��� ����� ������� ��� ����
    w_solved(i,1:size(w_solved_str)) =sort(double(w_solved_str));%��� �����, ������� ����������� -- ��������� ���� ���������
end

cr_val=1000000;%����� ����������� ������� ��� ��������, ���� � ��� ���� ������ ����� -- ����� ����� ��� ������ �� ���� ���������
for i=1:1:length(k)
    for j=1:1:size(w_solved_str)
        if imag(w_solved(i,j))~=0||(w_solved(i,j)==0)
            w_solved(i,j)=cr_val;
        end
    end
end

figure;
for i=1:1:length(k)
          j=1;
      if w_solved(i,j)~=cr_val %���� � ������� ��� ������ �����, �� �������� �����
           hold on; plot(k(i),w_solved(i,j),'.c','LineWidth',2);
      end
          j=2;
      if w_solved(i,j)~=cr_val %���� � ������� ��� ������ �����, �� �������� �����
           hold on; plot(k(i),w_solved(i,j),'.r','LineWidth',2);
      end
          j=3;
      if w_solved(i,j)~=cr_val %���� � ������� ��� ������ �����, �� �������� �����
           hold on; plot(k(i),w_solved(i,j),'.b','LineWidth',2);
      end
          j=4;
      if w_solved(i,j)~=cr_val %���� � ������� ��� ������ �����, �� �������� �����
           hold on; plot(k(i),w_solved(i,j),'.m','LineWidth',2);
      end
    
end
ylabel('$\omega$')
xlabel('$k$')

