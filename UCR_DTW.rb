=begin
  
  Author: Ire Sun
  Date: Dec 13, 2014
  Ref: Searching and Mining Trillions of Time Series Subsequences under Dynamic Time Warping, 2011

=end

class Trillion

  def self.dist (x, y)
    return (x-y)**2
  end

  def self.upper_lower_lemire (t, len, r)
    du = Array.new()
    dl = Array.new()
    u = Array.new(len)
    l = Array.new(len)

    du.push(0)
    dl.push(0)

    (1 .. len-1).each{ |i|

      if (i > r)
        u[i-r-1] = t[du.first.to_i]
        l[i-r-1] = t[dl.first.to_i]
      end

      # Pop out the bound that is not maximium or minimum.
      # Store the max upper bound and min lower bound within window r.
      if (t[i] > t[i-1])
        du.pop
        while (!(du.empty?) && t[i] > t[du.last.to_i])
          du.pop
        end
      else
        dl.pop
        while (!(dl.empty?) && t[i] < t[dl.last.to_i])
          dl.pop
        end
      end
      du.push(i)
      dl.push(i)

      # i - r - 1 == r + du.first
      # Pop out the bound that os out of window r.
      if ( i == 2 * r + 1 + du.first)
        du.shift
      elsif ( i == 2 * r + 1 + dl.first)
        dl.shift
      end
    }

    # The envelop of first r points are from r+1 .. r+r, so the last r points' envelop haven't settle down yet.
    (len .. len+r).each{ |i|
      u[i-r-1] = t[du.first.to_i]
      l[i-r-1] = t[dl.first.to_i]
      if (i - du.first.to_i >= 2 * r + 1)
        du.shift
      end
      if (i - dl.first.to_i >= 2 * r + 1)
        dl.shift
      end
    }

    return l, u
  end

  # O(1)
  def self.lb_kim_hierarchy(t, q, j, mean, std, bsf = Float::INFINITY)
    lb = 0

    # 1st points at front and back
    x0 = (t[j] - mean) / std
    y0 = (t[(j + q.length - 1)] - mean) / std
    lb = dist(x0, q.first) + dist(y0, q.last)
    return lb, 1 if( lb >= bsf )

    # 2 points at front
    x1 = (t[(j + 1)] - mean) / std
    d = [dist(x1, q.first), dist(x0, q[1]), dist(x1, q[1])].min
    lb += d
    return lb, 2 if ( d >= bsf)
    return lb, 1 if (lb >= bsf)

      # 2 points at back
      y1 = (t[(j + q.length - 2)] - mean) / std
      d = [dist(y1, q.last), dist(y0, q[-2]), dist(y1, q[-2])].min
      lb += d
    return lb, q.length - 2 + 1 if ( d >= bsf)
      return lb, 1 if (lb >= bsf)

      # 3 points at front
      x2 = (t[(j + 2)] - mean) / std
      d = [dist(x0,q[2]), dist(x1, q[2]), dist(x2,q[2]), dist(x2,q[1]), dist(x2,q.first)].min
      lb += d
    return lb, 3 if ( d >= bsf)
      return lb, 1 if (lb >= bsf)

      # 3 points at back
      y2 = (t[(q.length - 3 + j)] - mean) / std
      d = [dist(y0, q[-3]), dist(y1, q[-3]), dist(y2, q[-3]), dist(y2, q[-2]), dist(y2, q[-1])].min
      lb += d
    return lb, q.length - 3 + 1 if ( d >= bsf)

      return lb, 1
  end

  # O(n)
  def self.lb_keogh_cumulative(order, t, uo, lo, cb, j, m, mean, std, bsf = Float::INFINITY)
    lb = 0
    jump = order.first

    # Compute the distances of t and query envelop.
    m.times{ |i|
      x = (t[(j + order[i]).to_i] - mean) / std
      d = 0

      jump = order[i] if(order[i] < jump)

      if (x > uo[i])
        d = dist(x, uo[i])
      elsif (x < lo[i])
        d = dist(x, lo[i])
      end
      lb += d
      cb[order[i]] =  d

      if(lb >= bsf)
        break
      end
    }

    return lb, (jump + 1)
  end

  # O(n)
  def self.lb_keogh_data_cumulative (order, qo, cb, l, u, j, mean, std, bsf = Float::INFINITY)
    lb = 0
    jump = order.first

    order.length.times{ |i|
      uu = (u[j + order[i]] - mean) / std
      ll = (l[j + order[i]] - mean) / std
      d = 0

      jump = order[i] if(order[i] < jump)

      if (qo[i] > uu)
        d = dist(qo[i], uu)
      elsif (qo[i] < ll)
        d = dist(qo[i], ll)
      end
      lb += d
      cb[order[i]] = d

      if(lb >= bsf)
        break
      end
    }
    return lb, (jump + 1)
  end

  # Draw the diagonally dtw window path to a rectangle, so we can reuse the cost/cost_prev array.
  def self.dtw(a_seq, b_seq, cb, r, bsf = Float::INFINITY)
    cost = Array.new(2 * r + 1, Float::INFINITY)
    cost_prev = Array.new(2 * r + 1, Float::INFINITY)
    k = 0

    a_seq.length.times{ |i|
      k = [0, r - i].max
      min_cost = Float::INFINITY

      j = [0,i-r].max
      while(j <= [a_seq.length-1,i+r].min)
        if( i == 0 && j == 0)
          min_cost = cost[k] = dist(a_seq[0], b_seq[0])
        else
          if ((j-1<0)||(k-1<0))
            y = Float::INFINITY
                else
                  y = cost[k-1]
                end
                if ((i-1<0)||(k+1>2*r))
                  x = Float::INFINITY
                else
            x = cost_prev[k+1]
                end
                if ((i-1<0)||(j-1<0))
                  z = Float::INFINITY
                else
                  z = cost_prev[k]
                end

                cost[k] = [x, y, z].min + dist(a_seq[i], b_seq[j])

                if(cost[k] < min_cost)
                  min_cost = cost[k]
                end
        end
        k += 1
        j += 1
      end

      if(i+r < a_seq.length-1 && min_cost + cb[i+r+1] >= bsf)
        return min_cost + cb[i+r+1]
      end

      cost_tmp = cost
          cost = cost_prev
          cost_prev = cost_tmp
    }

    k -= 1

    return cost_prev[k]
  end

  EPOCH = 100000

  def calcuate(dat, query, window_rate, dontSort, dontJump)

    loc = 0
    kim = keogh = keogh2 = 0
    dist = lb_kim = lb_k = lb_k2 = 0.0

    data_name = dat
    query_name = query

    ex = ex2 = 0.0
    bsf = Float::INFINITY
    i = j = 0

    start = Time.now

    # Read query
    puts "Read query"
    qp = File.open(query_name, "r")
    q = Array.new()
    qp.read.split(" ").each{ |x| 
      q.push(x.to_f)
      ex += q.last
      ex2 += q.last**2
    }
    qp.close

    mean = ex / q.length
    std = (ex2 / q.length - mean ** 2) ** 0.5

    # Normalize
    puts "Normalize query"
    q.map!{ |x|
      (x - mean) / std
    }

    if(window_rate <= 1)  # wraping window
      r = (window_rate*q.length).floor
    else
      r = window_rate.floor
    end

    # Calculating the envelop of q
    puts "Calculating the envelop of q"
    l, u = Trillion::upper_lower_lemire(q, q.length, r)

    q_tmp = q.each_with_index.map{ |x, i|
      IndexObj.new(i, x)
    }

    order = Array.new()
    qo = Array.new()
    lo = Array.new()
    uo = Array.new()

    if(!dontSort)
      puts "Sorting"

      q_tmp.sort!{ |a, b|
        (b.value.abs) <=> (a.value.abs)
      }

      q_tmp.each{ |x|
        order.push(x.index)
        qo.push(q[x.index])
        uo.push(u[x.index])
        lo.push(l[x.index])
      }
    else
      order = (0 .. q.length - 1).to_a
      qo = q
      uo = u
      lo = l
    end

    cb = Array.new(q.length, 0)
    cb1 = Array.new(q.length, 0)
    cb2 = Array.new(q.length, 0)

    i = j  = jump_times = 0
    ex = ex2 = 0.0
    done = false
    it = ep = k = 0

    tp = File::open(data_name, "r")
    tdata = tp.read.split(" ").map{ |x|
      x.to_f
    }
    tp.close

    buffer = Array.new(EPOCH)
    t = Array.new(q.length * 2)
    tz = Array.new(q.length * 2)

    puts "Searching"

    while(!done)
      ep = 0
      if (it == 0)
        (q.length-1).times{ |k|
          if(tdata.length != 0)
            buffer[k] = tdata.shift
          end
        }
      else
        (q.length-1).times{ |k|
          buffer[k] = buffer[EPOCH - q.length + 1 + k]
        }
      end

      # Read a chunk(=EPOCH) into buffer
      ep = q.length - 1
      while(ep < EPOCH && tdata.length != 0)
        buffer[ep] = tdata.shift
        ep += 1
      end

      if(ep <= q.length - 1)
        done = true
      else
        l_buff, u_buff = Trillion::upper_lower_lemire(buffer, ep, r)

        if (it % (1000000/(EPOCH - q.length + 1)) == 0 )
          print "."
        end

        ex = ex2 = 0
        jump = 0
        ep.times{|i|
          d = buffer[i]

          ex += d
          ex2 += d ** 2

          t[i % q.length] = d
          t[(i % q.length) + q.length] = d

          jump -= 1 if( jump > 0 )

          # After the size of t is q.length or more, starts the pruning works
          if( i >= q.length - 1)
            j = (i + 1) % q.length

            if( dontJump || jump <= 0 )
              mean = ex / q.length
              std = (ex2 / q.length - mean ** 2) ** 0.5

              iCap = i - (q.length - 1)

              lb_kim, jump = Trillion::lb_kim_hierarchy(t, q, j, mean, std, bsf)

              if(lb_kim < bsf)
                lb_k, jump = Trillion::lb_keogh_cumulative(order, t, uo, lo, cb1, j, q.length, mean, std, bsf)

                if(lb_k < bsf)

                  lb_k2, jump = Trillion::lb_keogh_data_cumulative(order, qo, cb2, l_buff, u_buff, iCap, mean, std, bsf)

                  if (lb_k2 < bsf)

                    if(lb_k > lb_k2)
                      cb[q.length - 1] = cb1[q.length - 1]
                      (q.length - 2).downto(0){|k|
                        cb[k] = cb[k + 1] + cb1[k]
                      }
                    else
                      cb[q.length - 1] = cb2[q.length - 1]
                      (q.length - 2).downto(0){|k|
                        cb[k] = cb[k + 1] + cb2[k]
                      }
                    end

                    tz = t[j .. j + q.length - 1].map{ |x|
                      (x - mean) / std
                    }

                    dist = Trillion::dtw(tz, q, cb, r, bsf)

                    if(dist < bsf)
                      bsf = dist
                      loc = it * (EPOCH - q.length + 1) + i - q.length + 1
                    end

                  else
                    keogh2 += 1
                  end
                else
                  keogh += 1
                end
              else
                kim += 1
              end
            else
              jump_times += 1
            end
            ex -= t[j]
            ex2 -= t[j] ** 2
          end
        }

        if (ep < EPOCH)
          done = true
        else
          it += 1
        end
      end

    end

    finish = Time.now

    puts " "

    @location = loc
    puts "Location: #{loc}"
    dist = bsf ** 0.5
    puts "Distance: #{dist}"
    i = it * (EPOCH - q.length + 1) + ep
    puts "Data scanned: #{i}"
    exeTime = finish - start
    puts "Execution time: #{exeTime}"
    puts " "

    jump_perc = 100.0 * jump_times / i
    puts "Pruned by Jump: #{jump_perc}%"
    kim_perc = 100.0 * kim / i
    puts "Pruned by LB_Kim: #{kim_perc}%"
    keogh_perc = 100.0 * keogh / i
    puts "Pruned by LB_Keogh: #{keogh_perc}%"
    keogh2_perc = 100.0 *keogh2 / i
    puts "Pruned by LB_Keogh2: #{keogh2_perc}%"
    dtw_prec = 100.0 - (kim_perc + keogh_perc + keogh2_perc + jump_perc)
    puts "DTW Calcuation  : #{dtw_prec}%"
  end

  def location
    @location
  end
end

class IndexObj
  def initialize(i, val)
    @index = i
    @value = val
  end

  def index
    @index
  end

  def value
    @value
  end
end

if __FILE__ == $PROGRAM_NAME
  if (ARGV.length < 1 || ARGV[0] == nil)
    puts "Need data file name as 1st parameter."
    exit
  end

  if (ARGV.length < 2 || ARGV[1] == nil)
    puts "Need query file name as 2nd parameter."
    exit
  end

  if (ARGV.length < 3 || ARGV[2] == nil)
    puts "Need window rate as 3rd parameter."
    exit
  end

  dontSort = false
  dontJump = false

  if (ARGV.length > 3)
    ARGV[3 .. ARGV.length - 1].each{ |x|
      dontSort = true if(x == "-ns")
      dontJump = true if(x == "-nj")
    }
  end

  Trillion.new(ARGV[0], ARGV[1], ARGV[2].to_f, dontSort, dontJump)

end