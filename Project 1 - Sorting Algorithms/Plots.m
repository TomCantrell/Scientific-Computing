
% Plots for sorting algorithms

bubble_sort = importdata("Bubble_sort_time.txt");
quick_sort = importdata("Quick_sort_time.txt");
heap_sort = importdata("Heap_sort_time.txt");

% Log plot for Bubblesort
figure(1)
plot(bubble_sort(:,1),bubble_sort(:,2)*10^-3)
xlabel("n")
ylabel("sort time / seconds")
grid on
%title("Bubblesort Algorithm")

% Log plot for Quicksort
figure(2)
plot(quick_sort(:,1),quick_sort(:,2)*10^-3)
xlabel("n")
ylabel("sort time / seconds")
grid on
%title("Quicksort Algorithm")

% Log plot for Heapsort
figure(3)
plot(heap_sort(:,1),heap_sort(:,2)*10^-3)
xlabel("n")
ylabel("sort time / seconds")
grid on
%title("Heapsort Algorithm")

% All algorithms
figure(4)
plot(bubble_sort(:,1),log(bubble_sort(:,2)*10^-3),quick_sort(:,1),log(quick_sort(:,2)*10^-3),heap_sort(:,1),log(heap_sort(:,2)*10^-3))
title("All algorithms")
legend("Bubblesort","Quicksort","Heapsort","Location", "North West")
xlabel("n")
ylabel("log(sort time)")
grid on

% All algorithms
figure(5)
plot(bubble_sort(:,1),bubble_sort(:,2)*10^-3,quick_sort(:,1),quick_sort(:,2)*10^-3,heap_sort(:,1),heap_sort(:,2)*10^-3)
title("All algorithms")
legend("Bubblesort","Quicksort","Heapsort","Location", "North West")
xlabel("n")
ylabel("sort time /seconds")
grid on

quicksort_ratio = [];
heapsort_ratio = [];
nlogn = [];
for i=2:length(quick_sort(:,2))
    nlogn = [nlogn;2*log(2)/log(quick_sort(i,1)) + 2 ];
    quicksort_ratio = [quicksort_ratio;quick_sort(i,2)/quick_sort(i-1,2)];
    heapsort_ratio = [heapsort_ratio;heap_sort(i,2)/heap_sort(i-1,2)];
end
figure(6)
x = 1:length(nlogn);
plot(x,quicksort_ratio,"--",x,nlogn,"-o")
grid on
legend("Factor calculated from quicksort data","Theoretical factor")

ylabel("Factor of increase in av. sort time")
figure(7)
plot(x,heapsort_ratio,"--",x,nlogn,"-o")
grid on
legend("Factor calculated from heapsort data","Theoretical factor","Location","SouthEast")
ylabel("Factor of increase in av. sort time")

figure(8)
plot(x,heapsort_ratio,"--",x,quicksort_ratio,"-.",x,nlogn,"-o")
legend("factor calculated from heapsort data","factor calculated from quicksort data","theoretical factor","Location","SouthEast")
grid on
ylabel("Factor of increase in av. sort time")

%%% Testing randomness of data
v = importdata("test_random.txt");
s_n = [];
vect = v(:,1);
for i=1:length(vect)
    x = vect(i);
    sum=0;
    for j=1:length(vect)
       if vect(j) <= x
           sum = sum + 1;
       end
    end
    s_n = [s_n; sum/length(vect)];
end
plot(s_n,".")
grid on
xlabel("vector entry")
ylabel("Proportion vector entries <= x value")














